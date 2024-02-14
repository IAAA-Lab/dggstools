import shutil
import tempfile
import unittest
from pathlib import Path

import rasterio.enums
import rhealpixdggs.dggs as rhp
import time

from sqlalchemy import create_engine

from dggstools.auids.auids import AUIDGenerator
from dggstools.auids.rhpx_auids import optimize_cuids_rhealpix, is_optimal_cuids_rhealpix, RHEALPixAUIDGenerator
from dggstools.rhpx.rhpxdataframes import *
from dggstools.rhpx.vector_to_rhpx import _vector_to_optimal_set_of_cuids
from dggstools.rhpx.vector_to_rhpx import *
from dggstools.rhpx.utils import utils
from dggstools.rhpx.utils.storage import geodataframe_to_geopackage, get_gpkg_rhpx_metadata, rhealpix_to_geopackage, \
    geopackage_to_rhealpix


# TODO: Many tests here are essentially smoke tests and the results produced by them must be checked manually if you
# need to be sure that they produce the correct results. Curating a set of test datasets, with the expected results
# after the different transformations, and using that in the tests to make sure that the results are what is
# expected is pending work.
class DataTestsSpec(unittest.TestCase):

    # This runs once for the whole class
    @classmethod
    def setUpClass(cls):
        cls.rdggs = rhp.RHEALPixDGGS(ellipsoid=WGS84_ELLIPSOID, north_square=1, south_square=0, N_side=3)
        cls.rdggs_helper = RHEALPixDGGSHelper(cls.rdggs)
        cls.logger = logging.getLogger("geo2dggs")
        cls.rdggs_df_helper = RHEALPixDataFrameHelper(cls.rdggs)
        cls.engine = create_engine("sqlite://", echo=True, future=True)

        # If I wanted to not see, say, rasterio logs:
        # logging.getLogger('rasterio').propagate = False
        # Configuration for the logs (including those of every library used that produces logs)
        logging.basicConfig(format="[%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s", level=logging.INFO)

        cwd = os.getcwd()
        try:
            data_dir = os.environ['GEO2DGGS_TEST_DATA_DIR']
            cls.data_dir = data_dir
        except KeyError:
            if os.path.exists("../test_data"):
                cls.data_dir = Path("../test_data")
            elif os.path.exists("test_data"):
                cls.data_dir = Path("test_data")
            elif os.path.exists("tests/test_data"):
                cls.data_dir = Path("tests/test_data")
        os.chdir(cwd)
        cls.temp_dir = tempfile.mkdtemp()
        shutil.copytree(cls.data_dir, cls.temp_dir, dirs_exist_ok=True)
        cls.data_dir = Path(cls.temp_dir)
        cls.output_dir = cls.temp_dir
        os.chdir(cls.temp_dir)

    def log_some_info(self, input_file_path, dst_resolution_idx):
        with rasterio.open(input_file_path) as raster:
            # Take the input file top, left, bottom and right coordinates of the input raster
            left, top, right, bottom, res_x, res_y = get_bbox_from_raster_profile(raster.profile)
            self.logger.info(f"BBOX OF: {[left, top, right, bottom, res_x, res_y]}")
            self.logger.info(f"CRS: {raster.profile['crs']}")
            top_left_cell = self.rdggs.cell_from_point(dst_resolution_idx, (left, top))
            bottom_right_cell = self.rdggs.cell_from_point(dst_resolution_idx, (right, bottom))
            if top_left_cell is not None:
                self.logger.info(f"TOP LEFT CELL centroid: {[top_left_cell.suid, top_left_cell.centroid()]}")
                self.logger.info(f"TOP LEFT CELL upper left: {[top_left_cell.suid, top_left_cell.ul_vertex()]}")
            else:
                self.logger.info("TOP LEFT CELL is None")

            if bottom_right_cell is not None:
                self.logger.info(f"BOTTOM RIGHT CELL centroid: {[bottom_right_cell.suid, bottom_right_cell.centroid()]}")
                self.logger.info(f"BOTTOM RIGHT upper left: {[bottom_right_cell.suid, bottom_right_cell.ul_vertex()]}")
            else:
                self.logger.info("BOTTOM RIGHT CELL is None")

    def _test_raster(self, input_file_name, dst_resolution_idx, input_crs, nc_variable, resampling):
        output_file_name = utils.insert_suffix(input_file_name, "-RHEALPIX")
        output_file_name = utils.change_extension(output_file_name, "tif")

        input_file_path = os.path.join(self.data_dir, input_file_name)
        if nc_variable is not None:
            output_file_name = utils.insert_suffix(output_file_name, f"-{nc_variable}")
            input_file_path = "netcdf:" + input_file_path + f":{nc_variable}"

        print(input_file_path)
        res_idx = raster_to_rhealpix(self.rdggs, input_file_path, os.path.join(self.output_dir, output_file_name),
                                     dst_resolution_idx, input_crs=input_crs, resampling=resampling)
        self.assertEqual(res_idx, dst_resolution_idx)
        self.log_some_info(os.path.join(self.output_dir, output_file_name), dst_resolution_idx)


    def _test_vector(self, input_file_name, dst_resolution_idx, input_crs, property_name, layer):
        output_file_name = utils.change_extension(input_file_name, "tif")

        res_idx, props_to_ints = vector_to_rhealpix(self.rdggs, os.path.join(self.data_dir, input_file_name),
                                                    os.path.join(self.output_dir, output_file_name), dst_resolution_idx,
                                                    input_crs=input_crs,
                                                    property_for_class=property_name, layer=layer)
        self.assertEqual(res_idx, dst_resolution_idx)
        self.log_some_info(os.path.join(self.output_dir, output_file_name), dst_resolution_idx)

    def test_raster_landsat_small(self):
        input_file_name = "landsat_image_small.tif"
        output_file_name = utils.insert_suffix(input_file_name, "-RHEALPIX")
        input_file_path = os.path.join(self.data_dir, input_file_name)
        output_file_path = os.path.join(self.output_dir, output_file_name)
        output_file_path_from_cells = utils.change_extension(utils.insert_suffix(output_file_path, "-fromgdf"), "tif")

        # After rescaling to a higher level than the one necessary to not lose any information from input file
        # the rmse and bias should be zero
        res_idx = raster_to_rhealpix(self.rdggs, input_file_path, output_file_path,
                                     rescaling_strategy=RescalingStrategy.TO_HIGHER)
        rmse, bias = calculate_rmse(input_file_path, output_file_path, every_rows=5, every_cols=5,
                                                        band=1)
        self.assertEqual(rmse, 0)
        self.assertEqual(bias, 0)
        self.logger.info(f"Between {input_file_path} and  {output_file_path} RMSE: {rmse}, BIAS: {bias}")


        gdf = self.rdggs_df_helper.rhealpix_file_to_geodataframe(output_file_path)
        self.rdggs_df_helper.geodataframe_to_rhealpix_file(gdf, output_file_path_from_cells)
        # After generating the cell values for the rhealpix raster dataset and then using them to initialize a new raster
        # those should be equal. For now I will just check that the rmse, bias, and the inverses, are all 0
        rmse, bias = calculate_rmse(output_file_path, output_file_path_from_cells, every_rows=5, every_cols=5,
                                                        band=1)
        inv_rmse, inv_bias = calculate_rmse(output_file_path_from_cells, output_file_path, every_rows=5, every_cols=5,
                                                        band=1)
        self.assertEqual(rmse, 0)
        self.assertEqual(bias, 0)
        self.assertEqual(inv_rmse, 0)
        self.assertEqual(inv_bias, 0)
        self.logger.info(f"Between {output_file_path} and {output_file_path_from_cells} RMSE: {rmse}, BIAS: {bias}, INV RMSE: {inv_rmse}, INV BIAS: {inv_bias}")

    # Tests for now are pretty naive. Unless there are some exceptions, they will pass
    def test_raster_datasets(self):
        datasets = [
            # TODO: After reprojecting the nasadem file to rhealpix and showing in QGIS it has a weird behaviour:
            # for example some rows seem to have smaller (less "tall") pixels than others. There might be some issues
            # with QGIS visualization (visualization of rhealpix some times is tricky), but this needs to be checked.
            ("nasadem.tif", 9, None, None, rasterio.enums.Resampling.nearest),
            ("Spain-France.tif", 7, None, None, rasterio.enums.Resampling.nearest),
            ("pnoa_2015_25830_0354_4_4.jpg", 7, rasterio.crs.CRS().from_epsg(25830), None, rasterio.enums.Resampling.nearest),
            ("NAIP_30.img", 7, None, None, rasterio.enums.Resampling.nearest),
            ("HYP_50M_SR_Robinson.tif", 7, None, None, rasterio.enums.Resampling.nearest),
            ("HYP_50M_SR_W.tif", 7, None, None, rasterio.enums.Resampling.cubic),
            ("Iberia.tif", 7, None, None, rasterio.enums.Resampling.cubic),
            ("c_gls_SCE500_202112070000_CEURO_MODIS_V1.0.1.nc", 8, None, None, rasterio.enums.Resampling.nearest),
        ]

        for input_file_name, res_idx, input_crs, nc_variable, resampling in datasets:
            self._test_raster(input_file_name, res_idx, input_crs, nc_variable, resampling)

    def _test_raster_alignments(self, input_file_name, dst_resolution_idx, input_crs):
        output_file_name = utils.insert_suffix(input_file_name, "-RHEALPIX-ALIGNRHP")
        output_file_name = utils.change_extension(output_file_name, "tif")

        res_idx = raster_to_rhealpix(self.rdggs, os.path.join(self.data_dir, input_file_name),
                                     os.path.join(self.output_dir, output_file_name), dst_resolution_idx,
                                     input_crs=input_crs)

        self.log_some_info(os.path.join(self.output_dir, output_file_name), dst_resolution_idx)

        gdf = self.rdggs_df_helper.rhealpix_file_to_geodataframe(os.path.join(self.output_dir,
                                                            output_file_name), add_uid=True)
        geodataframe_to_geopackage(gdf, os.path.join(self.output_dir,
                                                     utils.change_extension(output_file_name, "gpkg")))
        self.rdggs_df_helper.geodataframe_to_rhealpix_file(gdf,
                                      os.path.join(self.output_dir,
                                                   utils.insert_suffix(output_file_name, "-FROMGEODATAFRAME")))

    def test_raster_alignments(self):
        datasets = [
            ("Spain-France.tif", 7, None),
            ("landsat_image_small.tif", 8, None),
            ("pnoa_2015_25830_0354_4_4.jpg", 7, rasterio.crs.CRS.from_epsg(25830)),
            ("NAIP_30.img", 7, None),
            ("HYP_50M_SR_W.tif", 4, None),
            ("c_gls_SCE500_202112070000_CEURO_MODIS_V1.0.1.nc", 6, None),
        ]

        for input_file_name, res_idx, input_crs in datasets:
            self._test_raster_alignments(input_file_name, res_idx, input_crs)

    # TODO: Support for some of these datasets axis combination is not yet ready in vectorutils. Fix that
    # def test_vector_datasets(self):
    #     datasets = [
    #         ("Comunidades_Autonomas_ETRS89_30N.shp", 7, None, "Codigo", None),
    #         ("Aragón_ETRS89_30N.shp", 8, rasterio.crs.CRS.from_epsg(25830), "CODIGO", None),
    #         ("Pop_AgricRegion.shp", 7, None, "OBJECTID", None),
    #         ("ne_110m_land.shp", 4, None, None, None),
    #         ("antarctica.gpkg", 7, None, None, None),
    #         ("oprvrs_gb_WGS84.gpkg", 9, None, None, None)
    #     ]
    #
    #     for input_file_name, res_idx, input_crs, property_name, layer in datasets:
    #         self._test_vector(input_file_name, res_idx, input_crs, property_name, layer)


    def test_calculate_rmse(self):
        raster_to_rhealpix(self.rdggs, os.path.join(self.data_dir, "landsat_image_small.tif"),
                           os.path.join(self.output_dir, "landsat_image_small-RHEALPIX-RMSE.tif"), 10,
                           resampling=rasterio.enums.Resampling.bilinear)
        rmse, bias = calculate_rmse(os.path.join(self.data_dir, "landsat_image_small.tif"),
                                    os.path.join(self.output_dir, "landsat_image_small-RHEALPIX-RMSE.tif"),
                                    every_rows=30, every_cols=30, band=1)
        self.logger.info(f"RMSE: {rmse}, BIAS: {bias} landsat_image_small.tif -> landsat_image_small-RHEALPIX-RMSE.tif")

        # Inverses (taking as "truth"/control the RHEALPIX version)
        rmse, bias = calculate_rmse(os.path.join(self.output_dir, "landsat_image_small-RHEALPIX-RMSE.tif"),
                                    os.path.join(self.data_dir, "landsat_image_small.tif"),
                                    every_rows=30, every_cols=30, band=1)
        self.logger.info(f"RMSE: {rmse}, BIAS: {bias} landsat_image_small-REHALPIX-RMSE.tif -> landsat_image_small.tif")


        # This tif originally does not have properly set the nodata value. The RMSE is catastrofically badly calculated if
        # it assumes that -32768 is a proper value, so the best option is setting the nodata in the file
        # ourselves
        with rasterio.open(os.path.join(self.data_dir, "nasadem.tif"), "r+") as nasadem:
            nasadem.nodata = -32768

        raster_to_rhealpix(self.rdggs, os.path.join(self.data_dir, "nasadem.tif"),
                           os.path.join(self.output_dir, "nasadem-RHEALPIX-RMSE.tif"), 9,
                           resampling=rasterio.enums.Resampling.nearest)
        rmse, bias = calculate_rmse(os.path.join(self.data_dir, "nasadem.tif"),
                                    os.path.join(self.output_dir, "nasadem-RHEALPIX-RMSE.tif"),
                                    every_rows=30, every_cols=30, band=1)
        self.logger.info(f"RMSE: {rmse}, BIAS: {bias} nasadem.tif -> nasadem-RHEALPIX-RMSE.tif")

        # Inverses (taking as "truth"/control the RHEALPIX version)
        rmse, bias = calculate_rmse(os.path.join(self.output_dir, "nasadem-RHEALPIX-RMSE.tif"),
                                    os.path.join(self.data_dir, "nasadem.tif"),
                                    every_rows=30, every_cols=30, band=1)
        self.logger.info(f"RMSE: {rmse}, BIAS: {bias} nasadem-RHEALPIX-RMSE.tif -> nasadem.tif")


    def _test_calculate_vector_raster_area_error(self, input_file_name, dst_resolution_idx, input_crs, property_name,
                                                 layer):
        input_path = os.path.join(self.data_dir, input_file_name)
        output_path = os.path.join(self.output_dir, utils.change_extension(
            utils.insert_suffix(input_file_name, "-AREAERROR"), "tif"))
        output_path_tmp_file = os.path.join(self.output_dir,
                                            utils.change_extension(
                                                utils.insert_suffix(input_file_name, "-AREAERROR_tmp"), "tif"))

        _, _ = vector_to_rhealpix(self.rdggs, input_path, output_path,
                                  input_crs=input_crs,
                                  dst_resolution_idx=dst_resolution_idx,
                                  property_for_class=property_name,
                                  layer=layer)

        error, bias, vec_area, cell_area = calculate_vector_raster_area_error(
            input_path, output_path, property_for_class=property_name, layer=layer)


        self.logger.info(
            f"******************* Error {input_file_name}: {error / 1000000} km², bias {bias / 1000000} km²")
        self.logger.info(
            f"******************* Areas {input_file_name}: Vec: {vec_area / 1000000} km², cells: {cell_area / 1000000} km²")
        self.logger.info(
            f"******************* Ratio Error/Vec Area {input_file_name}: {error / vec_area}")

    # TODO: Support for some of these datasets axis combination is not yet ready in vectorutils. Fix that
    # def test_calculate_vector_raster_area_error(self):
    #     datasets = [
    #         ("Aragón_ETRS89_30N.shp", 8, None, "CODIGO", None),
    #         ("Pop_AgricRegion.shp", 8, None, "OBJECTID", None),
    #         ("antarctica.gpkg", 7, None, "surface", None),
    #         ("Comunidades_Autonomas_ETRS89_30N.shp", 7, None, "Codigo", None),
    #         ("Pop_AgricRegion.shp", 7, None, "OBJECTID", None),
    #     ]
    #
    #     for input_file_name, dst_resolution_idx, input_crs, property_name, layer in datasets:
    #         self._test_calculate_vector_raster_area_error(input_file_name, dst_resolution_idx, input_crs, property_name,
    #                                                       layer)

    def _test_calculate_vector_raster_line_error(self, input_file_name, dst_resolution_idx, input_crs, property_name,
                                                 layer, every_feature):
        input_path = os.path.join(self.data_dir, input_file_name)
        output_path = os.path.join(self.output_dir, utils.change_extension(
            utils.insert_suffix(input_file_name, "-LINEERROR"), "tif"))
        output_path_tmp_file = os.path.join(self.output_dir,
                                            utils.change_extension(
                                                utils.insert_suffix(input_file_name, "-LINEERROR_tmp"), "tif"))

        # Lines are better rasterized with this function instead of vector_to_rhealpix
        _, _ = vector_to_rhealpix_without_intermediate_raster(self.rdggs, input_path, output_path,
                                                              input_crs=input_crs,
                                                              dst_resolution_idx=dst_resolution_idx,
                                                              property_for_class=property_name,
                                                              layer=layer)

        error_per_node, error_per_feature = calculate_vector_raster_line_error(
            input_path, output_path, property_for_class=property_name, every_feature=every_feature, layer=layer)

        self.logger.info(
            f"******************* Error {input_file_name}: mean distance per node {error_per_node / 1000} km, "
            f"mean distance per feature {error_per_feature / 1000} km.")

    # TODO: Support for some of these datasets axis combination is not yet ready in vectorutils. Fix that
    # def test_calculate_vector_raster_line_error(self):
    #     datasets = [
    #         ("oprvrs_gb_WGS84.gpkg", 9, None, None, None, 200),  # With res 12 the process is killed in my
    #         # computer (32 GB of RAM), I suppose that because it needs too much memory.
    #     ]
    #
    #     for input_file_name, dst_resolution_idx, input_crs, property_name, layer, every_feature in datasets:
    #         self._test_calculate_vector_raster_line_error(input_file_name, dst_resolution_idx, input_crs,
    #                                                       property_name, layer, every_feature)

    def test_rhealpix_to_gpkg(self):
        self._test_raster("landsat_image_small.tif", 8, None, None, rasterio.enums.Resampling.nearest)
        gdf = self.rdggs_df_helper.rhealpix_file_to_geodataframe(os.path.join(self.output_dir, "landsat_image_small-RHEALPIX.tif"))

        geodataframe_to_geopackage(gdf, os.path.join(self.output_dir, "landsat_image_small.gpkg"))
        gdf.to_parquet(os.path.join(self.output_dir, "landsat_image_small.parquet"))
        gdf.to_feather(os.path.join(self.output_dir, "landsat_image_small.feather"))
        self.rdggs_df_helper.geodataframe_to_rhealpix_file(gdf, os.path.join(self.output_dir, "landsat_image_small_fromgdf.tiff"))

        self._test_raster("c_gls_SCE500_202112070000_CEURO_MODIS_V1.0.1.nc", 6, None, None, rasterio.enums.Resampling.nearest)
        gdf = self.rdggs_df_helper.rhealpix_file_to_geodataframe(os.path.join(self.output_dir,
                                                         "c_gls_SCE500_202112070000_CEURO_MODIS_V1.0.1-RHEALPIX.tif"))

        geodataframe_to_geopackage(gdf, os.path.join(self.output_dir,
                                                     "c_gls_SCE500_202112070000_CEURO_MODIS_V1.0.1.gpkg"))
        self.rdggs_df_helper.geodataframe_to_rhealpix_file(gdf, os.path.join(self.output_dir,
                                                   "c_gls_SCE500_202112070000_CEURO_MODIS_V1.0.1_fromgdf.tiff"))

        self._test_raster("HYP_50M_SR_Robinson.tif", 4, None, None, rasterio.enums.Resampling.nearest)
        gdf = self.rdggs_df_helper.rhealpix_file_to_geodataframe(os.path.join(self.output_dir,
                                                         "HYP_50M_SR_Robinson-RHEALPIX.tif"),
                                            add_uid=True)

        geodataframe_to_geopackage(gdf, os.path.join(self.output_dir,
                                                     "HYP_50M_SR_Robinson.gpkg"))
        self.rdggs_df_helper.geodataframe_to_rhealpix_file(gdf,
                                      os.path.join(self.output_dir,
                                                   "HYP_50M_SR_Robinson_fromgdf.tiff"))

    def test_rhealpix_to_gdf_with_json(self):
        self._test_raster("landsat_image_small.tif", 8, None, None, rasterio.enums.Resampling.nearest)
        gdf = self.rdggs_df_helper.rhealpix_file_to_geodataframe(
            os.path.join(self.output_dir, "landsat_image_small-RHEALPIX.tif"), values_in_json=True)
        geodataframe_to_geopackage(gdf, os.path.join(self.output_dir,
                                                     "landsat_image_small-json-RHEALPIX.gpkg"))

    def test_rhealpix_to_gdf_update_with_json(self):
        self._test_raster("landsat_image_small.tif", 8, None, None, rasterio.enums.Resampling.nearest)
        output_file_path = os.path.join(self.output_dir, "landsat_image_small-RHEALPIX.tif")
        gdf = self.rdggs_df_helper.rhealpix_file_to_geodataframe(output_file_path, values_in_json=True,
                                                                 store_nodata=True)
        geodataframe_to_geopackage(gdf, os.path.join(self.output_dir,
                                                     "landsat_image_small-json-RHEALPIX.gpkg"))
        # Modify the rhealpix raster file
        with rasterio.open(output_file_path, "r") as raster:
            new_band_1 = np.zeros(shape=(raster.width, raster.height))
            kwargs = raster.meta.copy()
            output_file_path_1 = os.path.join(self.output_dir, utils.insert_suffix("landsat_image_small-RHEALPIX.tif", "-1"))
            with rasterio.open(output_file_path_1, 'w', **kwargs) as dst:
                dst.write(new_band_1, 1)

        start_time_ns = time.process_time_ns()
        gdf1 = self.rdggs_df_helper.rhealpix_file_to_geodataframe_update(output_file_path_1, gdf, values_in_json=True)
        end_time_ns = time.process_time_ns()
        logger.info(f"Update geodataframe is {(end_time_ns - start_time_ns) / 1e6} ms")

        start_time_ns = time.process_time_ns()
        gdf2 = self.rdggs_df_helper.rhealpix_file_to_geodataframe(output_file_path_1, values_in_json=True)
        end_time_ns = time.process_time_ns()
        logger.info(f"Re-create geodataframe is {(end_time_ns - start_time_ns) / 1e6} ms")

        geodataframe_to_geopackage(gdf1, os.path.join(self.output_dir,
                                                     "landsat_image_small-json-RHEALPIX-1.gpkg"))
        geodataframe_to_geopackage(gdf2, os.path.join(self.output_dir,
                                                      "landsat_image_small-json-RHEALPIX-2.gpkg"))

    def test_gpkg_metadata(self):
        # Generate gpkg
        self._test_raster("landsat_image_small.tif", 9, None, None,
                          rasterio.enums.Resampling.nearest)
        output_file_path = os.path.join(self.output_dir, "landsat_image_small-RHEALPIX.tif")
        gdf = self.rdggs_df_helper.rhealpix_file_to_geodataframe(output_file_path, values_in_json=True,
                                                                 store_nodata=True)
        geodataframe_to_geopackage(gdf, os.path.join(self.output_dir,
                                                     "landsat_image_small-json-metadata-RHEALPIX.gpkg"))
        metadata = get_gpkg_rhpx_metadata(os.path.join(self.output_dir,
                                                     "landsat_image_small-json-metadata-RHEALPIX.gpkg"))
        self.assertEqual(9, metadata["res_idx"])
        self.assertEqual(8, metadata["nbands"])



    def test_generate_rhealpix_grid(self):
        # I want cells of approximately 100.000 meters of width/height
        # idx, res = self.rdggs_helper.get_closest_resolution(100000)
        # Or I want a fixed resolution index
        idx=0
        # TOFIX: This does not work for idx=5 and other resolutions, with different bounding boxes.
        # This appears to be due to some vertices from some cells falling out of the rhealpix "image".
        # I am not sure where the problem is: some vertices of some cells in rHEALPix are actually
        # outside of those cells. But the vertices operation in the Cell class, should take/is taking
        # that into account? Is other thing?
        # A "patch" is avoiding those checks. This is done in the next lines

        import rhealpixdggs.pj_healpix, rhealpixdggs.pj_rhealpix

        good_in_healpix_image = rhealpixdggs.pj_healpix.in_healpix_image
        good_in_rhealpix_image = rhealpixdggs.pj_rhealpix.in_rhealpix_image

        # We temporarily "patch" the rhealpix library to make sure that all cells are
        # included in the grid. This is not the best solution, but allows to generate all grids
        # (although it is not clear if every cell in them is correct, or if some cells
        # are missing. TODO: Find out if there is a definitive best solution for this
        rhealpixdggs.pj_healpix.in_healpix_image = lambda x, y: True
        rhealpixdggs.pj_rhealpix.in_rhealpix_image = lambda x, y, south_square, north_square: True

        gdf = self.rdggs_df_helper.rhealpix_grid_as_geodataframe((-180, 90), (180, -90), idx, as_geodetic=True)
        # As GeoPackage
        output_file_path = os.path.join(self.output_dir, "rhealpix_cells")
        geodataframe_to_geopackage(gdf, output_file_path + ".gpkg")
        # As GeoJSON
        gdf.to_file(output_file_path + ".json", driver="GeoJSON")

        rhealpixdggs.pj_healpix.in_healpix_image = good_in_healpix_image
        rhealpixdggs.pj_rhealpix.in_rhealpix_image = good_in_rhealpix_image

    def test_auid_generation_for_shape(self):
        input_file_name = "Aragón_ETRS89_30N.shp"
        input_crs = None
        dst_resolution_idx = 6
        property_name = "CODIGO"
        layer=None

        input_path = os.path.join(self.data_dir, input_file_name)
        output_path = os.path.join(self.output_dir, utils.change_extension(
            utils.insert_suffix(input_file_name, "-AUID"), "tif"))

        _, _ = vector_to_rhealpix(self.rdggs, input_path, output_path,
                                  input_crs=input_crs,
                                  dst_resolution_idx=dst_resolution_idx,
                                  property_for_class=property_name,
                                  layer=layer)

        gdf = self.rdggs_df_helper.rhealpix_file_to_geodataframe(output_path, add_uid=True)
        # Keep just Aragon shape: band1 == 2 (0 is the background)
        gdf = gdf[gdf['band1'] == 2]

        cuids =gdf['cellid'].values

        # Measuring the execution time in process time (not including the time spent in the OS)
        start_time_ns = time.process_time_ns()
        a = AUIDGenerator()
        auid_comp_b64, _ = a.generate_auid_hash_b64(cuids)
        cuids = a.cuids_from_auid_b64(auid_comp_b64)
        end_time_ns = time.process_time_ns()
        assert (sorted(cuids) == cuids)  # As long as there are no repetitions in cuids this must be true
        logger.info(f"AUID generation time for Aragón at resolution {dst_resolution_idx} "
                    f"is {(end_time_ns - start_time_ns)/1e6} ms")

    def test_auid_optimizer_for_shape(self):
        input_file_name = "Aragón_ETRS89_30N.shp"
        input_crs = None
        dst_resolution_idx = 8
        property_name = "CODIGO"
        layer = None

        input_path = os.path.join(self.data_dir, input_file_name)
        output_path = os.path.join(self.output_dir, utils.change_extension(
            utils.insert_suffix(input_file_name, "-AUID"), "tif"))

        _, _ = vector_to_rhealpix(self.rdggs, input_path, output_path,
                                  input_crs=input_crs,
                                  dst_resolution_idx=dst_resolution_idx,
                                  property_for_class=property_name,
                                  layer=layer)

        gdf = self.rdggs_df_helper.rhealpix_file_to_geodataframe(output_path, add_uid=True)
        # Keep just Aragon shape: band1 == 2 (0 is the background)
        gdf = gdf[gdf['band1'] == 2]

        cuids = gdf['cellid'].values

        a = AUIDGenerator()
        auid_comp_b64, _ = a.generate_auid_hash_b64(cuids)
        cuids = a.cuids_from_auid_b64(auid_comp_b64)
        # Measuring the execution time in process time (not including the time spent in the OS)
        start_time_ns = time.process_time_ns()
        optimal_cuids = optimize_cuids_rhealpix(cuids, n_side=self.rdggs.N_side)
        end_time_ns = time.process_time_ns()
        auid_optimal_comp_b64, _ = a.generate_auid_hash_b64(optimal_cuids)
        logger.info(f"AUID optimization time for Aragón at resolution {dst_resolution_idx} "
                    f"is {(end_time_ns - start_time_ns) / 1e6} ms")
        logger.info(f"Optimal AUID is {100*len(auid_optimal_comp_b64)/len(auid_comp_b64)}% of the size of the base "
                    f"and has {100*len(optimal_cuids)/len(cuids)}% of the number of cuids of the base.")

    def test_optimal_auid_generation_for_shape(self):

        # Execution times (quick tests, just for reference) for Aragón:
        # Res 6
        # Optimal creation: 2s
        # Optimized a posteriori: 0,5 s
        # Res 7
        # Optimal creation: 13s
        # Optimized a posteriori: 4s
        # Res 8
        # Optimal creation: 120s
        # Optimized a posteriori: 33s ("slow" algorithm was 72s)
        # Res 9
        # Optimal creation: 1200s
        # Optimized a posteriori: 349 s ("slow" algorithm was **5000s**)
        #
        # So, the times seem to grow exponentially with the resolution (approx 9~10 times longer for each resolution
        # which is reasonable because n_side 3, so 9 times more cells per resolution increment). Optimized a posteriori
        # is much faster for large numbers of cells.

        input_file_name = "Aragón_ETRS89_30N.shp"
        input_crs = None
        dst_resolution_idx = 7
        property_name = "CODIGO"
        layer = None

        output_base_file_path = os.path.join(self.output_dir,
                                             change_extension(insert_suffix(input_file_name, "-AUID"), "tif"))

        # Measuring the execution time in process time (not including the time spent in the OS)
        start_time_ns = time.process_time_ns()
        optimal_cuids = _vector_to_optimal_set_of_cuids(self.rdggs, os.path.join(self.data_dir, input_file_name),
                                                        output_base_file_path, input_crs=input_crs,
                                                        dst_resolution_idx=dst_resolution_idx,
                                                        property_for_class=property_name,
                                                        included_classes=[2], layer=layer)

        end_time_ns = time.process_time_ns()
        logger.info(f"Optimal generation: {sorted(optimal_cuids)}")
        logger.info(f"Optimal generation is {(end_time_ns - start_time_ns)/1e6} ms in creating optimal.")

        # To be fair with times, I should rasterize again the highest resolution, but that is not really
        # taking too long, so the times are still useful
        start_time_ns = time.process_time_ns()
        highest_res_path = insert_suffix(output_base_file_path, f"-{dst_resolution_idx}")

        # Let's compare with the highest resolution raster optimized a posteriori
        gdf = self.rdggs_df_helper.rhealpix_file_to_geodataframe(highest_res_path, add_uid=True)
        # Keep just Aragon shape: band1 == 2 (0 is the background)
        gdf = gdf[gdf['band1'] == 2]
        cuids = gdf['cellid'].values

        geodataframe_to_geopackage(gdf, highest_res_path + ".gpkg")

        a = AUIDGenerator()
        auid_comp_b64, _ = a.generate_auid_hash_b64(cuids)
        cuids = a.cuids_from_auid_b64(auid_comp_b64)
        optimized_cuids = optimize_cuids_rhealpix(cuids, n_side=self.rdggs.N_side)
        end_time_ns = time.process_time_ns()

        logger.info(f"Optimized a posteriori cuids: {sorted(optimized_cuids)}")
        logger.info(f"Optimización a posteriori, partiendo del rasterizado highest res, is {(end_time_ns - start_time_ns) / 1e6} ms")

        assert(sorted(optimized_cuids) == sorted(optimal_cuids))

        start_time_ns = time.process_time_ns()
        assert(is_optimal_cuids_rhealpix(optimal_cuids, n_side=self.rdggs.N_side))
        end_time_ns = time.process_time_ns()
        logger.info(f"Optimality test on the result takes: {(end_time_ns - start_time_ns) / 1e6} ms")

        start_time_ns = time.process_time_ns()
        optimize_cuids_rhealpix(optimized_cuids, n_side=self.rdggs.N_side)
        end_time_ns = time.process_time_ns()
        logger.info(f"Optimizing when it is not necessary takes: {(end_time_ns - start_time_ns) / 1e6} ms")

    def test_generate_mini_rhealpix_from_df(self):
        # We generate some "uncommon" cell combinations to test that they are correctly handled
        # Some generated files may later be used to test other things as they are small and thus very fast
        output_file_path = os.path.join(self.temp_dir, "mini_rhealpix_n1_s0_nside3_Ncap.tiff")

        df = DataFrame([{"entity_id": "N0", "value": 0}])
        df = df.astype({"value": np.int16})

        gdf = self.rdggs_df_helper.rhealpix_data_frame_to_geodataframe(df, nodata=-9999, value_columns=["value"],
                                                                       geo_id_column_name="entity_id")
        self.rdggs_df_helper.geodataframe_to_rhealpix_file(gdf, output_file_path, nodata=-9999)
        logger.info(f"Generated {output_file_path}.")
        with rasterio.open(output_file_path) as raster:
            self.assertEqual(raster.nodata, -9999)
            self.assertEqual(raster.width, 1)
            self.assertEqual(raster.height, 1)
            self.assertEqual(raster.count, 1)

        df = DataFrame([{"entity_id": "N0", "value": 0},
                        {"entity_id": "N1", "value": 1},
                        {"entity_id": "N2", "value": 2},
                        {"entity_id": "N3", "value": 3},
                        {"entity_id": "N4", "value": 4},
                        {"entity_id": "N5", "value": 5},
                        {"entity_id": "N6", "value": 6},
                        {"entity_id": "N7", "value": 7},
                        {"entity_id": "N8", "value": 8},
                        ])
        df = df.astype({"value": np.int16})

        gdf = self.rdggs_df_helper.rhealpix_data_frame_to_geodataframe(df, nodata=-9999, value_columns=["value"],
                                                                       geo_id_column_name="entity_id")
        self.rdggs_df_helper.geodataframe_to_rhealpix_file(gdf, output_file_path, nodata=-9999)
        logger.info(f"Generated {output_file_path}.")
        with rasterio.open(output_file_path) as raster:
            self.assertEqual(raster.nodata, -9999)
            self.assertEqual(raster.width, 3)
            self.assertEqual(raster.height, 3)
            self.assertEqual(raster.count, 1)

        output_file_path = os.path.join(self.temp_dir, "mini_rhealpix_n1_s0_nside3_Scap.tiff")
        df = DataFrame([{"entity_id": "S0", "value": 0},
                        {"entity_id": "S1", "value": 1},
                        {"entity_id": "S2", "value": 2},
                        {"entity_id": "S3", "value": 3},
                        {"entity_id": "S4", "value": 4},
                        {"entity_id": "S5", "value": 5},
                        {"entity_id": "S6", "value": 6},
                        {"entity_id": "S7", "value": 7},
                        {"entity_id": "S8", "value": 8},

                        ])
        df = df.astype({"value": np.int16})

        gdf = self.rdggs_df_helper.rhealpix_data_frame_to_geodataframe(df, nodata=-9999, value_columns=["value"],
                                                                       geo_id_column_name="entity_id")
        self.rdggs_df_helper.geodataframe_to_rhealpix_file(gdf, output_file_path, nodata=-9999)
        logger.info(f"Generated {output_file_path}.")
        with rasterio.open(output_file_path) as raster:
            self.assertEqual(raster.nodata, -9999)
            self.assertEqual(raster.width, 3)
            self.assertEqual(raster.height, 3)
            self.assertEqual(raster.count, 1)

        output_file_path = os.path.join(self.temp_dir, "mini_rhealpix_n1_s0_nside3_N_and_S.tiff")
        df = DataFrame([{"entity_id": "S0", "value": 0},
                        {"entity_id": "S1", "value": 1},
                        {"entity_id": "S2", "value": 2},
                        {"entity_id": "S3", "value": 3},
                        {"entity_id": "S4", "value": 4},
                        {"entity_id": "S5", "value": 5},
                        {"entity_id": "S6", "value": 6},
                        {"entity_id": "S7", "value": 7},
                        {"entity_id": "S8", "value": 8},
                        {"entity_id": "N0", "value": 0},
                        {"entity_id": "N1", "value": 1},
                        {"entity_id": "N2", "value": 2},
                        {"entity_id": "N3", "value": 3},
                        {"entity_id": "N4", "value": 4},
                        {"entity_id": "N5", "value": 5},
                        {"entity_id": "N6", "value": 6},
                        {"entity_id": "N7", "value": 7},
                        {"entity_id": "N8", "value": 8},
                        ])
        df = df.astype({"value": np.int16})

        gdf = self.rdggs_df_helper.rhealpix_data_frame_to_geodataframe(df, nodata=-9999, value_columns=["value"],
                                                                       geo_id_column_name="entity_id")
        self.rdggs_df_helper.geodataframe_to_rhealpix_file(gdf, output_file_path, nodata=-9999)
        logger.info(f"Generated {output_file_path}.")
        with rasterio.open(output_file_path) as raster:
            self.assertEqual(raster.nodata, -9999)
            # 6 because the rdggs parameters n1 and s0 (+ possibly 1 col of nodata cells accounting for some rounding)
            self.assertTrue(raster.width == 6 or raster.width == 7)
            # 9 because the rdggs ; + possibly 1 row of nodata cells accounting for some rounding)
            self.assertTrue(raster.height == 10 or raster.height == 9)
            self.assertEqual(raster.height, 10) #
            self.assertEqual(raster.count, 1)

    def test_raster_to_gpkg_and_back(self):
        # Reproject a raster to rhpx so we have something to work with
        _ = raster_to_rhealpix(self.rdggs,  os.path.join(self.data_dir, "landsat_image_small.tif"),
                               os.path.join(self.temp_dir, "landsat_image_small_rhpx_321.tif"),
                               dst_resolution_idx=9)

        rhealpix_to_geopackage(os.path.join(self.temp_dir, "landsat_image_small_rhpx_321.tif"),
                               os.path.join(self.temp_dir, "landsat_image_small_tif_321.gpkg"))

        geopackage_to_rhealpix(os.path.join(self.temp_dir, "landsat_image_small_tif_321.gpkg"),
                               os.path.join(self.temp_dir, "landsat_image_small_gpkg_321.tif"))

        # At this point landsat_image_small_rhpx_321.tif and landsat_image_small_gpkg_321.tif should be equivalent
        # save from minor details (e.g. compression scheme in the GeoTIFF). TODO: TEST THIS, for now this have
        # been manually verified, but it should be automated.



if __name__ == '__main__':
    unittest.main()
