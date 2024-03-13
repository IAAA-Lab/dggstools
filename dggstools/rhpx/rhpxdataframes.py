import json

import numpy
import sys
from pandas import DataFrame

from .raster_to_rhpx import *
from .rhpxutils import *
from .utils.rasterutils import *
from .utils.vectorutils import bounds_to_left_top_right_bottom
from geopandas import GeoDataFrame
from .utils.utils import almost_equal
import uuid

logger = logging.getLogger(__name__)


class _NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.integer):
            return int(obj)
        if isinstance(obj, numpy.floating):
            return float(obj)
        if isinstance(obj, numpy.ndarray):
            return obj.tolist()
        if isinstance(obj, uuid.UUID):
            return str(obj)
        return super(_NpEncoder, self).default(obj)

class RHEALPixDataFrameHelper:
    """
    """
    def __init__(self, rdggs: RHEALPixDGGS):
        self.rdggs = rdggs
        self.rdggs_helper = RHEALPixDGGSHelper(self.rdggs)


    def rhealpix_file_to_geodataframe(self, input_file_path: str,
                                      geo_id_column_name: str = "cellid", add_uid: bool = False,
                                      values_in_json: bool = False, store_nodata: bool = False) -> GeoDataFrame:
        """
        Takes a raster file with the CRS corresponding to self.rdggs, and generates a GeoDataFrame with the contents
        of its non empty cells, explicit cellids and some metadata.
        """
        profile = get_raster_profile(input_file_path)
        left, top, right, bottom, resx, resy = get_bbox_from_raster_profile(profile)

        resolution_idx_x, _ = self.rdggs_helper.get_closest_resolution(abs(resx))
        resolution_idx_y, _ = self.rdggs_helper.get_closest_resolution(abs(resy))  # resy is often a negative number
        assert resolution_idx_x == resolution_idx_y, \
            f"{input_file_path} is not a proper rhealpix file. Its cells are not squares."

        resolution_idx = resolution_idx_x

        # It seems that cells_from_region fails if the rectangle defined by the two points passed includes
        # "empty space". I.e., if it includes a part outside any side of the cube. For example a rectangle
        # that includes part of the N or S squares, and part of the O,P,Q or R rectangles, and has a part
        # which is outside.
        # So this does not work in the general case:
        # return rdggs.cells_from_region(resolution_idx, top_left_cell.ul_vertex(), bottom_right_cell.ul_vertex())

        # Generate the cell ids (and values after that) contained in a given rhealpix file without
        # using the cells_from_region (or improve that function)

        cells = {geo_id_column_name: [], "geometry": []}
        if values_in_json:
            cells["all_bands"] = []
        if add_uid:
            cells["uuid"] = []

        with rasterio.open(input_file_path) as raster:

            logger.info(f"raster.nodata {raster.nodata}")
            logger.info(f"raster.nodatavals {raster.nodatavals}")

            if not values_in_json:
                for nband in range(raster.count):
                    cells[f"band{nband+1}"] = []   # We label band columns starting in 1 (as in gdal)
            data = raster.read()
            is_aligned = True
            all_pixels_transformed = True
            for row in range(raster.height):
                for col in range(raster.width):
                    # Read all the bands in row,col
                    cell_data = [data[i, row, col] for i in range(raster.count)]
                    # If any band in a cell is different from nodata, we need to store that cell
                    # If all were nodata, we wouldn't have to store that cell unless store_nodata is True
                    if any([cd != nd for cd, nd in zip(cell_data, raster.nodatavals)]) or store_nodata:
                        x, y = raster.xy(row, col)
                        # Getting the cell for each point is slow, but safer than using cells_from_region in the
                        # general case
                        curr_cell = self.rdggs.cell_from_point(resolution_idx, (x, y))
                        if curr_cell is not None:
                            cells[geo_id_column_name].append(str(curr_cell))
                            if not values_in_json:
                                for nband in range(raster.count):
                                    cells[f"band{nband+1}"].append(cell_data[nband])
                            else:
                                all_bands_dict = {}
                                for nband in range(raster.count):
                                    all_bands_dict[f"band{nband+1}"] = cell_data[nband]
                                # If we don't use our _NpEncoder, types such as int16 are not properly converted to json
                                cells["all_bands"].append(json.dumps(all_bands_dict, cls=_NpEncoder))

                            curr_cell_centroid = curr_cell.centroid()
                            is_aligned = is_aligned and almost_equal(curr_cell_centroid[0], x) \
                                         and almost_equal(curr_cell_centroid[1], y)
                            cells['geometry'].append(shapely.geometry.Point(x, y))
                            if add_uid:
                                cells['uuid'].append(str(uuid.uuid4()))
                        else:
                            # This pixel must be outside the DGGS. If it had nodata value, we wouldn't be here
                            # but sometimes an RGB value of (0,0,0) is used instead of a proper nodata, and that
                            # can't be known for certain here. We will hope for the best and ignore this pixel
                            all_pixels_transformed = False
            if not is_aligned:
                logger.warning(f"{input_file_path} does not seem to be a properly aligned rhealpix raster. We are "
                               f"storing its x,y values with each cell as produced by rasterio in the geodataframe, "
                               f"but those do not correspond with cell centroids. You should check the output of this "
                               f"function, as it probably contains some repeated cellids (that would not be possible if "
                               f"every pixel in the input corresponded to exactly one cell).")
            if not all_pixels_transformed:
                logger.warning(f"Some pixels in {input_file_path} do not correspond to any cell in the DGGS. "
                               f"This is OK as long as those pixels have nodata values (you would not be reading this "
                               f"warning in that case), or values which should be nodata, but they are not, maybe "
                               f"because limitations in the file format (e.g., they could be RGB(0,0,0). You should check "
                               f"your dataset if you find problems with this transformation.")

        # The cellid column could be the index of the dataframe. Think about it.
        gdf = GeoDataFrame(cells, crs=self.rdggs_helper.rhealpixdef_to_proj_string())

        if not values_in_json:
            # This makes sure that the GeoDataFrame has the correct data types for the values. If we don't do this
            # some data types are not properly recognized (e.g., a uint8 may be considered an int8 etc.).
            conversions = {}
            for nband in range(raster.count):
                if raster.dtypes[nband] is not None:
                    conversions[f"band{nband+1}"] = raster.dtypes[nband]
            gdf = gdf.astype(conversions)

        gdf = self._modfiy_gdf_attrs(gdf, left, top, right, bottom, raster.height, raster.width, resx, resolution_idx,
                                        raster.count, raster.nodata, raster.nodatavals, store_nodata, raster.dtypes)

        return gdf

    def _modfiy_gdf_attrs(self, gdf: GeoDataFrame, left: float, top: float, right: float, bottom: float, height: int,
                          width: int, res: float, resolution_idx: int, nbands: int, nodata: float,
                          nodatavals: List[float], store_nodata: bool, dtypes: List[str]) -> GeoDataFrame:

        # attrs attribute is experimental at this time, be careful with its use as it may change. It might not be
        # automatically persisted when saving to geopackage, parquet etc.
        gdf.attrs["left"] = left
        gdf.attrs["top"] = top
        gdf.attrs["right"] = right
        gdf.attrs["bottom"] = bottom
        gdf.attrs["height"] = height
        gdf.attrs["width"] = width
        gdf.attrs["res"] = res
        gdf.attrs["res_idx"] = resolution_idx
        gdf.attrs["nbands"] = nbands
        # Even if we do not store the cells with nodata, it is interesting to store which values were
        # used for nodata in the original raster dataset, in case we are recreating it afterwards from this gdf.
        # Besides this, if some cell has data for some bands but not for others, we need to store something for that
        # cell in the bands with nodata
        gdf.attrs["nodata"] = nodata
        gdf.attrs["nodatavals"] = nodatavals
        gdf.attrs["store_nodata"] = store_nodata
        # Not every input format will have this, but it is worth storing it
        gdf.attrs["dtypes"] = dtypes
        # The CRS is stored along the GeoDataFrame, but some details of the parameters of the rhealpixdggs used are not
        # and we might be needing them
        ellipsoid = {}
        for (k, v) in sorted(self.rdggs.ellipsoid.__dict__.items()):
            if k == "phi_0":
                continue
            ellipsoid[k] = v
        gdf.attrs["rhealpixdggs"] = {"n_side": self.rdggs.N_side,
                                     "north_square": self.rdggs.north_square,
                                     "south_square": self.rdggs.south_square,
                                     "max_areal_resolution": self.rdggs.max_areal_resolution,
                                     "max_resolution": self.rdggs.max_resolution,
                                     "ellipsoid": ellipsoid}
        return gdf


    def rhealpix_file_to_geodataframe_update(self, input_file_path: str, original_gdf: GeoDataFrame,
                                             geo_id_column_name: str = "cellid",
                                             values_in_json: bool = False, store_nodata: bool = False) -> GeoDataFrame:
        """
        This is *much* faster than rhealpix_file_to_geodataframe. Useful if we are modifying a rhealpix raster and
        need to have it as a GeoDataFrame from time to time.
        """
        profile = get_raster_profile(input_file_path)
        left, top, right, bottom, resx, resy = get_bbox_from_raster_profile(profile)

        assert_msg = "You are trying to update a GeoDataFrame with an incompatible rhealpix raster " \
                     "(extent and/or resolution are different)."
        assert original_gdf.attrs["left"] == left, assert_msg
        assert original_gdf.attrs["top"] == top, assert_msg
        assert original_gdf.attrs["right"] == right, assert_msg
        assert original_gdf.attrs["bottom"] == bottom, assert_msg
        assert abs(original_gdf.attrs["res"]) == abs(resx), assert_msg
        assert abs(original_gdf.attrs["res"]) == abs(resy), assert_msg
        assert original_gdf.attrs["store_nodata"], "store_nodata must be True for the original GeoDataFrame or this " \
                                                   "method will not work in the general case."

        cells = original_gdf.to_dict()

        if values_in_json:
            if "all_bands" not in cells.keys():
                raise ValueError(f"original_gdf does not have a column named 'all_bands', but requires one if you set values_in_json to True")

        with rasterio.open(input_file_path) as raster:
            logger.info(f"raster.nodata {raster.nodata}")
            logger.info(f"raster.nodatavals {raster.nodatavals}")

            if not values_in_json:
                for nband in range(raster.count):
                    if f"band{nband + 1}" not in cells.keys():
                        raise ValueError(
                            f"original_gdf does not have a column named {nband + 1}, but requires one if you set values_in_json to False")
            data = raster.read()
            for row in range(raster.height):
                for col in range(raster.width):
                    # Read all the bands in row,col
                    cell_data = [data[i, row, col] for i in range(raster.count)]
                    # If any band in a cell is different from nodata, we need to store that cell
                    # If all were nodata, we wouldn't have to store that cell, unless store_nodata is True
                    if any([cd != nd for cd, nd in zip(cell_data, raster.nodatavals)]) or store_nodata:
                        idx = row * raster.width + col
                        if not values_in_json:
                            for nband in range(raster.count):
                                cells[f"band{nband + 1}"][idx] = cell_data[nband]
                        else:
                            all_bands_dict = {}
                            for nband in range(raster.count):
                                all_bands_dict[f"band{nband + 1}"] = cell_data[nband]
                            # If we don't use our _NpEncoder, types such as int16 are not properly converted to json
                            cells["all_bands"][idx] = json.dumps(all_bands_dict, cls=_NpEncoder)

        # The cellid column could be the index of the dataframe. Think about it.
        gdf = GeoDataFrame(cells, crs=self.rdggs_helper.rhealpixdef_to_proj_string())

        if not values_in_json:
            # This makes sure that the GeoDataFrame has the correct data types for the values. If we don't do this
            # some data types are not properly recognized (e.g., a uint8 may be considered an int8 etc.).
            for nband in range(raster.count):
                if raster.dtypes[nband] is not None:
                    gdf[f"band{nband + 1}"] = gdf[f"band{nband + 1}"].astype(raster.dtypes[nband])

        # attrs attribute is experimental at this time, be careful with its use as it may change. It might not be
        # automatically persisted when saving to geopackage, parquet etc.
        gdf.attrs = original_gdf.attrs
        gdf.attrs["store_nodata"] = store_nodata  # This can change from the original
        return gdf


    def geodataframe_to_rhealpix_file(self, gdf: GeoDataFrame, output_file_path: str,
                                      metadata_dict: dict = None,
                                      nodata: int | float = 0):
        """
        Takes a GeoDataFrame produced by the rhealpix_file_to_geodataframe method and produces a GeoTiff with
        the same contents.
        """
        output_crs = rasterio.crs.CRS.from_string(self.rdggs_helper.rhealpixdef_to_proj_string())

        logger.info(f"Input file bands, rows, cols and resolution_idx "
                    f"{[gdf.attrs['nbands'], gdf.attrs['height'], gdf.attrs['width'], gdf.attrs['res_idx']]}")
        logger.info(f"Input file left, right, top, bottom and resolution "
                    f"{[gdf.attrs['left']], gdf.attrs['right'], gdf.attrs['top'], gdf.attrs['bottom'], gdf.attrs['res']}")

        left = gdf.attrs["left"]
        right = gdf.attrs["right"]
        top = gdf.attrs["top"]
        bottom = gdf.attrs["bottom"]
        height = gdf.attrs["height"]
        width = gdf.attrs["width"]
        nbands = gdf.attrs["nbands"]
        res = gdf.attrs['res']

        # Create the right transform, given the rdggs and the data
        transform, transform_width, transform_height = rasterio.warp.calculate_default_transform(
            output_crs, output_crs, width, height,
            left=left, right=right, top=top, bottom=bottom)

        if transform[0] != res or transform[4] != -res:
            # Sometimes I have to force the output resolution or it will not be exact.
            # Forcing this, may change the width and/or height, so it is not the default choice;
            # it should not be a big issue: there will be some additional cells with nodata, but that's all.
            # It is not perfect, but as long as it passes the tests...
            transform, transform_width, transform_height = rasterio.warp.calculate_default_transform(
                output_crs, output_crs, width, height,
                left=left, right=right, top=top, bottom=bottom, resolution=(res,res))

        # If this is not true, the transform has not been as exact as expected
        assert(almost_equal(transform[0], res))

        dtype = gdf.attrs["dtypes"][0]  # First band dtype. If there are bands
        # with different dtypes, geotiff is not the right format (theoretical support, but not in practice)

        if dtype is None:
            raise ValueError("Without a given dtype in the geodataframe attrs, we can't generate the geotiff")
        if dtype not in ["byte", "uint8", "int8", "uint16", "int16", "uint32", "int32", "float32", "float64"]:
            raise ValueError(f"Unsupported dtype {dtype} for GeoTIFF files. "
                             f"See <https://gdal.org/drivers/raster/gtiff.html>.")

        # Create a numpy array of the right shape filled with the nodata value
        data_array = numpy.full((nbands, transform_height, transform_width), nodata, dtype=dtype)

        for tuple in gdf.itertuples():
            point = getattr(tuple, "geometry")
            row, col = rasterio.transform.rowcol(transform, point.x, point.y)
            for band in range(nbands):
                data_array[band, row, col] = getattr(tuple, f"band{band + 1}")

        # Generate output file
        kwargs = {
            'driver': 'GTiff',
            'compress': 'DEFLATE',
            'crs': output_crs,
            'transform': transform,
            'width': transform_width,
            'height': transform_height,
            'count': nbands,
            'dtype': dtype,
            'nodata': nodata
        }

        with rasterio.open(output_file_path, 'w', **kwargs) as dst:
            # To store the metadata in GeoTIFF, we use a tag (we name it user_metadata)
            dst.update_tags(user_metadata=json.dumps(metadata_dict, cls=_NpEncoder))
            for i in range(1, nbands+1):
                dst.write_band(i, data_array[i-1])

    def rhealpix_grid_as_geodataframe(self, nw: Tuple[int, int], se: Tuple[int, int], target_resolution: int,
                                      geo_id_column_name: str = "cellid",
                                      as_geodetic: bool = True) -> GeoDataFrame:
        """
        Generate all rhealpix cells:
         - With self.rdggs.
         - In the bounding box (ellipsoidal quadrangle to be more precise) defined by the north-west
         and south-east (nw, se) corners (both as (lon, lat) in geodetic coordinates).
         - At the given target_resolution (integer, min resolution is 0)
         - If as_geodetic, in EPSG:4326, if not, as rHEalPix (as defined by rdggs)
         - As a GeoDataFrame with two columns: entity_id_column_name and geometry. Geometries will typically be polygons with 4 corners, but
         if as_geodetic, they hay have less than 4, or they might even be linestrings (e.g., the N cap corners in geodetic
         coordinates are aligned).
        """
        cells = self.rdggs.cells_from_region(target_resolution, nw, se, plane=False)

        cell_data = {geo_id_column_name: [], "geometry": []}
        if as_geodetic:
            for row in cells:
                for cell in row:
                    cell_data[geo_id_column_name].append(str(cell))
                    cell_data['geometry'].append(shapely.geometry.Polygon(cell.vertices(plane=False,
                                                                                        trim_dart=True)))
            return GeoDataFrame(cell_data, crs="EPSG:4326")
        else:
            for row in cells:
                for cell in row:
                    cell_data[geo_id_column_name].append(str(cell))
                    cell_data['geometry'].append(shapely.geometry.Polygon(cell.vertices(plane=True)))
            return GeoDataFrame(cell_data, crs=self.rdggs_helper.rhealpixdef_to_proj_string())

    def rhealpix_data_frame_to_geodataframe(self, df: DataFrame, nodata: float = None,
                                            value_columns: List[str] = ["value"],
                                            geo_id_column_name: str = "cellid") -> GeoDataFrame:
        """
        PRE:
        - The df has a number of cells and all of them are at the same resolution level of the DGGS.
        - The df has one or more columns with values (n bands).
        - All of those columns have the same dtype and the same nodata (given as parameter).

        In general, this method produces a GeoDataFrame with the same structure (columns, metadata (attrs) etc.)
        as the one produced by rhealpix_file_to_geodataframe.

        Adds a geometry column to the df with the centroids of the cells which identifier is in geo_id_column_name
        and returns that as GeoDataFrame. Alsa adds the attrs in gdf_attrs (same format as produced by the
        get_gdf_attrs_from_rhealpix_file method) to the GeoDataFrame that can be deduced form the df, or
        given as parameter. And it will rename the value columns as "band1", "band2" etc.
        """

        nbands = len(value_columns)

        maxx = -sys.float_info.max
        maxy = -sys.float_info.max
        minx = sys.float_info.max
        miny = sys.float_info.max
        res_idx = None
        geometry = []
        for i in range(nbands):
            df.rename(columns={value_columns[i]: f"band{i+1}"}, inplace=True)  # band1 is the 0th element in the array

        for cellid in df[geo_id_column_name]:
            suid = self.rdggs.cell(cellidstr_to_suid(cellid))
            if res_idx is None:
                res_idx = cellid_resolution_idx(cellid)
            elif res_idx != cellid_resolution_idx(cellid):
                raise ValueError("All cells in the dataframe must have the same resolution")
            cell_minx, cell_miny, cell_maxx, cell_maxy = shapely.geometry.Polygon(suid.vertices(plane=True)).bounds
            minx = min(minx, cell_minx)
            miny = min(miny, cell_miny)
            maxx = max(maxx, cell_maxx)
            maxy = max(maxy, cell_maxy)
            geometry.append(shapely.geometry.Point(suid.centroid(plane=True)))

        left, top, right, bottom = bounds_to_left_top_right_bottom((minx, miny, maxx, maxy),
                                                                   self.rdggs_helper.rhealpixdef_to_proj_string())

        resx = self.rdggs.cell_width(res_idx)  # in the plane; in the ellipsoid the widths at the same resolution
                                               # level are different

        height = math.floor((maxy - miny) / resx) - 1 # resx must be equal to resy
        width = math.floor((maxx - minx) / resx) - 1

        # Prevents a division by zero. It is possible that this should be an exception as I don't imagine
        # a situation where this would be OK, but rounding and very low resolutions could lead to this, I suppose
        if height <= 0:
            height = 1
        if width <= 0:
            width = 1

        dtype = df["band1"].dtype  # We assume that the first column has the same dtype as the rest
        store_nodata = nodata is not None

        gdf = GeoDataFrame(df, crs=self.rdggs_helper.rhealpixdef_to_proj_string(), geometry=geometry)

        # Make sure that the calculated bounds are properly aligned to the dggs grid
        transform = self.rdggs_helper.align_bounds((left, top, right, bottom), res_idx)
        aligned_left = transform[2]
        aligned_top = transform[5]
        # The alignment should be "perfect" as it is based on the vertices of the cells
        assert(almost_equal(aligned_left, left))
        assert(almost_equal(aligned_top, top))

        # nodatavals = [nodata] -> PRE is that all nodata values in all the columns are the same
        gdf = self._modfiy_gdf_attrs(gdf, left=left, top=top, right=right, bottom=bottom,
                                     height=height, width=width, res=resx, resolution_idx=res_idx,
                                     nbands=nbands, nodata=nodata, nodatavals=[nodata], store_nodata=store_nodata,
                                     dtypes=[dtype])
        return gdf

