import unittest

import boto3
import rasterio.enums
import rhealpixdggs.dggs as rhp
from rasterio.session import AWSSession
from rhealpixdggs.ellipsoids import Ellipsoid

from dggstools.rhpx.rhpxdataframes import *
from dggstools.rhpx.vector_to_rhpx import *
from dggstools.rhpx.utils import utils
from dggstools.rhpx.utils.storage import geodataframe_to_geopackage, geodataframe_to_postgis


# TODO: Some proper tests of data, where the actual results are compared to a set of correct ones.
class ManualDataTestsSpec(unittest.TestCase):

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

    def _test_raster(self, input_file_name, dst_resolution_idx, input_crs, nc_variable=None,
                     resampling=rasterio.enums.Resampling.nearest):
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


    def setUp(self):
        self.data_dir = "/home/rbejar/Nextcloud_UNIZAR/Test_Data/input_files"
        self.output_dir = "/tmp"
        self.rdggs = rhp.RHEALPixDGGS(ellipsoid=WGS84_ELLIPSOID, north_square=1, south_square=0, N_side=3)
        self.rdggs_helper = RHEALPixDGGSHelper(self.rdggs)
        self.logger = logging.getLogger("geo2dggs")

        # If I wanted to not see, say, rasterio logs:
        # logging.getLogger('rasterio').propagate = False

        # Configuration for the logs (including those of every library used that produces logs)
        logging.basicConfig(format="[%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s", level=logging.INFO)

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

        rdggs_df_helper = RHEALPixDataFrameHelper(self.rdggs)

        gdf = rdggs_df_helper.rhealpix_file_to_geodataframe(output_file_path)
        rdggs_df_helper.geodataframe_to_rhealpix_file(gdf, output_file_path_from_cells)
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
        # We can change some parameters of the rdggs if we want

        self.rdggs = rhp.RHEALPixDGGS(ellipsoid=Ellipsoid(a=WGS84_A, f=WGS84_F, lon_0=10), north_square=1, south_square=0, N_side=3)
        datasets = [
            #("nasadem.tif", 9, None, None, rasterio.enums.Resampling.nearest),
            # ("Spain-France.tif", 7, None, None, rasterio.enums.Resampling.nearest),
            # ("landsat_image.tif", 7, None, None, rasterio.enums.Resampling.nearest),
            # ("pnoa_2015_25830_0354_4_4.jpg", 7, rasterio.crs.CRS.from_epsg(25830), None, rasterio.enums.Resampling.nearest),
            # ("NAIP_30.img", 7, None, None, rasterio.enums.Resampling.nearest),
            #("HYP_50M_SR.tif", 8, None, None, rasterio.enums.Resampling.nearest),
            #("HYP_50M_SR_W.tif", 7, None, None, rasterio.enums.Resampling.cubic),
            ("Iberia.tif", 9, None, None, rasterio.enums.Resampling.cubic),
            # ("c_gls_SCE500_202112070000_CEURO_MODIS_V1.0.1.nc", 9, None, None, rasterio.enums.Resampling.nearest),
            # ("c_gls_SWE5K_202204240000_NHEMI_SSMIS_V1.0.2.nc", 9, None, "swe", rasterio.enums.Resampling.nearest),
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

        output_file_name2 = utils.insert_suffix(input_file_name, "-RHEALPIX-ALIGNRASTERIO")
        output_file_name2 = utils.change_extension(output_file_name2, "tif")
        res_idx = raster_to_rhealpix(self.rdggs, os.path.join(self.data_dir, input_file_name),
                                     os.path.join(self.output_dir, output_file_name2), dst_resolution_idx,
                                     input_crs=input_crs)
        rdggs_df_helper = RHEALPixDataFrameHelper(self.rdggs)

        gdf = rdggs_df_helper.rhealpix_file_to_geodataframe(os.path.join(self.output_dir,
                                                         output_file_name), add_uid=True)
        geodataframe_to_geopackage(gdf, os.path.join(self.output_dir,
                                                     utils.change_extension(output_file_name, "gpkg")))
        rdggs_df_helper.geodataframe_to_rhealpix_file(gdf,
                                      os.path.join(self.output_dir,
                                                   utils.insert_suffix(output_file_name, "-FROMGEODATAFRAME")))

        gdf = rdggs_df_helper.rhealpix_file_to_geodataframe(os.path.join(self.output_dir,
                                                         output_file_name2), add_uid=True)
        geodataframe_to_geopackage(gdf, os.path.join(self.output_dir,
                                                     utils.change_extension(output_file_name2, "gpkg")))
        rdggs_df_helper.geodataframe_to_rhealpix_file(gdf, os.path.join(self.output_dir,
                                                   utils.insert_suffix(output_file_name2, "-FROMGEODATAFRAME")))


    def test_raster_alignments(self):
        datasets = [
            #("Spain-France.tif", 7, None),
            ("landsat_image.tif", 8, None),
            #("pnoa_2015_25830_0354_4_4.jpg", 7, rasterio.crs.CRS.from_epsg(25830)),
            #("NAIP_30.img", 7, None),
            #("HYP_50M_SR.tif", 5, None),
            #("c_gls_SCE500_202112070000_CEURO_MODIS_V1.0.1.nc", 9, None),
        ]

        for input_file_name, res_idx, input_crs in datasets:
            self._test_raster_alignments(input_file_name, res_idx, input_crs)

    def test_vector_datasets(self):
        # We can change lon_0 if we need to
        #self.rdggs = rhp.RHEALPixDGGS(ellipsoid=Ellipsoid(a=WGS84_A, f=WGS84_F, lon_0=10), north_square=1, south_square=0, N_side=3)
        datasets = [
            #("Comunidades_Autonomas_ETRS89_30N.shp", 7, None, "Codigo", None),
            ("Aragón_ETRS89_30N.shp", 9, rasterio.crs.CRS.from_epsg(25830), "CODIGO", None),
            #("Pop_AgricRegion.shp", 7, None, "OBJECTID", None),
            #("ne_110m_land.shp", 4, None, None, None),
            #("antarctica.gpkg", 7, None, None, None),
            #("oprvrs_gb.gpkg", 10, None, None, "WatercourseLink")
        ]

        for input_file_name, res_idx, input_crs, property_name, layer in datasets:
            self._test_vector(input_file_name, res_idx, input_crs, property_name, layer)

    def test_calculate_rmse(self):
        raster_to_rhealpix(self.rdggs, os.path.join(self.data_dir, "landsat_image.tif"),
                           os.path.join(self.output_dir, "landsat_image-RHEALPIX-RMSE.tif"), 10,
                           resampling=rasterio.enums.Resampling.bilinear)
        rmse, bias = calculate_rmse(os.path.join(self.data_dir, "landsat_image.tif"),
                                    os.path.join(self.output_dir, "landsat_image-RHEALPIX-RMSE.tif"),
                                    every_rows=30, every_cols=30, band=1)
        self.logger.info(f"RMSE: {rmse}, BIAS: {bias} landsat_image.tif -> landsat_image-RHEALPIX-RMSE.tif")

        # Inverses (taking as "truth"/control the RHEALPIX version)
        rmse, bias = calculate_rmse(os.path.join(self.output_dir, "landsat_image-RHEALPIX-RMSE.tif"),
                                    os.path.join(self.data_dir, "landsat_image.tif"),
                                    every_rows=30, every_cols=30, band=1)
        self.logger.info(f"RMSE: {rmse}, BIAS: {bias} landsat_image-REHALPIX-RMSE.tif -> landsat_image.tif")


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
    #
    #     # The errors are smaller (both visually, loading the datasets in QGIS, and measuring them) if I take
    #     # the original vector file (oprvrs_gb.gpkg, in EPSG:27700), transform it to EPSG:4326 and then transform
    #     # that one to rHEALPix instead of transforming directly the first one. I suppose that this could be
    #     # avoided with a more careful use of Proj, but for now I just will take notice. I don't seem to have
    #     # installed the best grid for the 27700 transformations in the Proj in my system, which is a bit old,
    #     # but I am not sure if the pyproj I am using here has it, or not). TODO: This issue is not a problem right
    #     # now, but it should be solved/fixed/documented or whatever appropriately
    #     # The transformation to WGS84 has been done with this function:
    #     # reproject_vector_file(os.path.join(self.data_dir, "oprvrs_gb.gpkg"),
    #     #                       os.path.join(self.output_dir, "oprvrs_gb_WGS84.gpkg"),
    #     #                       pyproj.CRS.from_epsg(4326), layer="WatercourseLink")
    #
    #     datasets = [
    #         ("oprvrs_gb_WGS84.gpkg", 9, None, None, "WatercourseLink", 200),  # With res 12 the process is killed in my
    #         # computer (32 GB of RAM), I suppose that because it needs too much memory.
    #         #("oprvrs_gb.gpkg", 11, None, None, "WatercourseLink", 50)
    #     ]
    #
    #     for input_file_name, dst_resolution_idx, input_crs, property_name, layer, every_feature in datasets:
    #         self._test_calculate_vector_raster_line_error(input_file_name, dst_resolution_idx, input_crs,
    #                                                       property_name, layer, every_feature)

    def test_rhealpix_to_gpkg(self):
        self._test_raster("landsat_image.tif", 9, None, None, rasterio.enums.Resampling.nearest)
        rdggs_df_helper = RHEALPixDataFrameHelper(self.rdggs)

        gdf = rdggs_df_helper.rhealpix_file_to_geodataframe(os.path.join(self.output_dir, "landsat_image-RHEALPIX.tif"))

        geodataframe_to_geopackage(gdf, os.path.join(self.output_dir, "landsat_image.gpkg"))
        gdf.to_parquet(os.path.join(self.output_dir, "landsat_image.parquet"))
        gdf.to_feather(os.path.join(self.output_dir, "landsat_image.feather"))
        rdggs_df_helper.geodataframe_to_rhealpix_file(gdf, os.path.join(self.output_dir, "landsat_image_fromgdf.tiff"))
        # landsat_image-REHALPIX.tif, resolution 9: Takes 7.5 seconds.
        # 214 KB in geotiff (zipping barely changes its size)
        # 465,4 KB parquet (if we zip the result, it shrinks to 278,1 KB)
        # 561,6 KB feather (if we zip the result, it shrinks to 264,2 KB)
        # 2,4 MB geopackage (if we zip the result, it shrinks to 507,4 KB). This is the only one that for now
        # includes metadata (gdf.attrs).

        # landsat_image-REHALPIX.tif, resolution 10: Takes 62 seconds.
        # 1,7 MB in geotiff
        # 3,6 MB parquet
        # 4,9 MB feather
        # 24,2 MB geopackage
        # includes metadata (gdf.attrs).

        # landsat_image-REHALPIX.tif, resolution 11: Takes 552 seconds (~9 minutes).
        # 14,5 MB in geotiff (zipping barely changes its size)
        # 31,2 MB parquet (if we zip the result, it shrinks to 278,1 KB)
        # 43,1 MB feather (if we zip the result, it shrinks to 264,2 KB)
        # 231,2 MB geopackage
        # includes metadata (gdf.attrs).

        # Processing times make sense. Each additional resolution level has 9 times more cells (with nside=3)
        # and the processing times _show that progression so far. Memory will grow the same, so resolutions of 12
        # and above (approx. 15 meters of resolution with that nside=3) will be difficult to process, for an
        # image this size, without chunking the data.
        # The times include: geotiff to gdf, gdf to gpkg, to parquet and to feather, and then gpkg to geotiff again

        # For this very small test with default compressions etc., it looks like parquet files are a good choice,
        # and zipping them is even better if trying to save bandwidth or storage space). GeoTiff is smaller of
        # course, but it lacks explicit cellids. Geopackage is a good format for testing, as it can be loaded
        # directly in QGIS to see the results. Feather might be still less stable than parquet and size is similar.
        # GEOPANDAS SUPPORT FOR PARQUET AND FEATHER AT THIS MOMENT IS NOT READY FOR PRODUCTION, AS THEY EXPECT
        # CHANGES IN THE GEO ASPECTS OF THOSE FORMATS. SEE THE LOGS THAT THOSE FUNCTIONS PRODUCE

        self._test_raster("c_gls_SCE500_202112070000_CEURO_MODIS_V1.0.1.nc", 6, None)
        gdf = rdggs_df_helper.rhealpix_file_to_geodataframe(os.path.join(self.output_dir,
                                                         "c_gls_SCE500_202112070000_CEURO_MODIS_V1.0.1-RHEALPIX.tif"))

        geodataframe_to_geopackage(gdf, os.path.join(self.output_dir,
                                                     "c_gls_SCE500_202112070000_CEURO_MODIS_V1.0.1.gpkg"))
        rdggs_df_helper.geodataframe_to_rhealpix_file(gdf, os.path.join(self.output_dir,
                                                   "c_gls_SCE500_202112070000_CEURO_MODIS_V1.0.1_fromgdf.tiff"))
        # c_gls_SCE500_202112070000_CEURO_MODIS_V1.0.1-REHALPIX.tif, resolution 9 (that is the closest higher resolution
        # to the original, it's ~50% better): Takes 4072 seconds (68 mins) (and some 7,5 GB of RAM in my laptop)
        # 3,9 MB in geotiff (the original .nc file is 172,4 MB, although it compresses to some 3 MB in ZIP)
        # 1,2 GB geopackage (231,MB zipped)
        # includes metadata (gdf.attrs).

        self._test_raster("HYP_50M_SR_Robinson.tif", 4, None)
        gdf = rdggs_df_helper.rhealpix_file_to_geodataframe(os.path.join(self.output_dir,
                                                         "HYP_50M_SR_Robinson-RHEALPIX.tif"),
                                            add_uid=True)

        geodataframe_to_geopackage(gdf, os.path.join(self.output_dir,
                                                     "HYP_50M_SR_Robinson.gpkg"))
        rdggs_df_helper.geodataframe_to_rhealpix_file(gdf,
                                      os.path.join(self.output_dir,
                                                   "HYP_50M_SR_Robinson_fromgdf.tiff"))

    # Only works with the proper s3 route and credentials
    def test_s3_support(self):
        input_file_name = "Spain-France.tif"
        input_file_path = "s3://data-tests/" + input_file_name
        output_file_name = utils.insert_suffix(input_file_name, "-RHEALPIX")
        output_file_path = "s3://data-tests/" + output_file_name

        # AWS S3 credentials must use boto3

        # In this case we are assuming that there is a ~/.aws/credentials file with a profile
        # named rids3 and with the proper aws_access_key_id and aws_secret_access_key values
        # (see <https://boto3.amazonaws.com/v1/documentation/api/latest/guide/credentials.html>)
        boto3_session = boto3.Session(profile_name='rids3')

        # This would also work, but it would hardcode the credentials here. It is OK for local testing though:
        # boto3_session = boto3.Session(aws_access_key_id='WHATEVER',
        #                              aws_secret_access_key='WH

        # I have my test aws_s3_endpoint in the same credentials file, but I have to read it as
        # any other ini file to retrieve it
        import configparser, pathlib
        cfg = configparser.ConfigParser()
        cfg.read(str(pathlib.Path.home())+'/.aws/credentials')
        aws_s3_endpoint = cfg['rids3']['aws_s3_endpoint']

        # Besides the AWSSession, other GDAL parameters are needed for S3 access to work, so we pass
        # them in a rasterio.Env. These parameters may be different for different S3 providers.
        # The GDAL documentation on virtual file systems (in particular vsis3) is a good starting point:
        # <https://gdal.org/user/virtual_file_systems.html>
        # CPL_VSIL_USE_TEMP_FILE_FOR_RANDOM_WRITE='YES' is required for writing GeoTIFF files to S3.
        # The AWS_S3_ENDPOINT is necessary for different S3 providers. For instance, for a local setup with
        # minio it could be something like "localhost:9000".
        with rasterio.Env(AWSSession(boto3_session),
                          AWS_S3_ENDPOINT=aws_s3_endpoint,
                          CPL_VSIL_USE_TEMP_FILE_FOR_RANDOM_WRITE='YES'):
            raster_to_rhealpix(self.rdggs, input_file_path, output_file_path)

    def test_generate_rhealpix_grid_gpkg(self):
        # We can change some parameters of the rdggs if we want
        self.rdggs = rhp.RHEALPixDGGS(ellipsoid=Ellipsoid(a=WGS84_A, f=WGS84_F, lon_0=10), north_square=1,
                                      south_square=0, N_side=3)
        rdggs_df_helper = RHEALPixDataFrameHelper(self.rdggs)
        # I want cells of approximately 100.000 meters of width/height
        #idx, res = get_closest_resolution(self.rdggs, 100000)
        #gdf = rhealpix_grid_as_geodataframe(self.rdggs, (-8, 48), (8, 32), idx, as_geodetic=True)
        gdf = rdggs_df_helper.rhealpix_grid_as_geodataframe((-180, 90), (180, -90), 1, as_geodetic=False)
        # As GeoPackage
        output_file_path = os.path.join(self.output_dir, "rhealpix_cells")
        geodataframe_to_geopackage(gdf, output_file_path + ".gpkg")
        # As GeoJSON
        #gdf.to_file(output_file_path + ".json", driver="GeoJSON")

    # Will just work if there is a PostgreSQL server running
    def test_generate_rhealpix_grid_postgre(self):
        # We can change some parameters of the rdggs if we want
        self.rdggs = rhp.RHEALPixDGGS(ellipsoid=Ellipsoid(a=WGS84_A, f=WGS84_F, lon_0=10), north_square=1,
                                         south_square=0, N_side=3)
        rdggs_df_helper = RHEALPixDataFrameHelper(self.rdggs)
        # I want cells of approximately 100.000 meters of width/height
        # idx, res = get_closest_resolution(self.rdggs, 100000)
        # gdf = rhealpix_grid_as_geodataframe(self.rdggs, (-8, 48), (8, 32), idx, as_geodetic=True)
        gdf = rdggs_df_helper.rhealpix_grid_as_geodataframe((-180, 90), (180, -90), 1, as_geodetic=False)

        # Only works with a PostgreSQL server running
        # def test_rhealpix_to_postgis(self):
        self._test_raster("landsat_image.tif", 9, None, None, rasterio.enums.Resampling.nearest)

        rdggs_df_helper = RHEALPixDataFrameHelper(self.rdggs)
        gdf = rdggs_df_helper.rhealpix_file_to_geodataframe(os.path.join(self.output_dir, "landsat_image-RHEALPIX.tif"))
        # It fails unless I do this (rhealpix_file_to_geodataframe creates a non-EPSG CRS, which does not seem supported
        # by sqlalchemy, which is behind geopandas.to_postgis)
        gdf.crs = None
        # I am testing locally with a Docker PostGIS with default values. Adapt it to your configuration
        # until we address properly the testing of this project
        geodataframe_to_postgis(gdf, "landsat_table", "mypg", "mypgpass", "localhost",
                                25432, "mypgdb", if_exists = "replace")

    # This fails, but I don't see why. Commenting out for now.
    def test_large_nc(self):
        input_file_name = "c_gls_SCE_202210250000_NHEMI_VIIRS_V1.0.1.nc"
        output_file_name = utils.insert_suffix(input_file_name, "-RHEALPIX")
        output_file_name = utils.change_extension(output_file_name, "tif")

        input_file_path = os.path.join(self.data_dir, input_file_name)

        res_idx = raster_to_rhealpix(self.rdggs, input_file_path, os.path.join(self.output_dir, output_file_name))
        self.log_some_info(os.path.join(self.output_dir, output_file_name), res_idx)



if __name__ == '__main__':
    unittest.main()
