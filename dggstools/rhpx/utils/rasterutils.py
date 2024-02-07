import math
import logging
from typing import Tuple, Union

import rasterio
import rasterio.crs
import rasterio.warp
import pyproj

logger = logging.getLogger(__name__)


def get_bbox_from_raster_profile(profile: rasterio.profiles.Profile):
    """
    Returns left, top, right, bottom res_x and res_y.
    res_y will normally (always?) be negative
    :param profile:
    :return:
    """
    left = profile["transform"][2]
    top = profile["transform"][5]
    res_x = profile["transform"][0]
    res_y = profile["transform"][4]
    right = left + profile["width"] * res_x
    bottom = top + profile["height"] * res_y
    return left, top, right, bottom, res_x, res_y


def get_raster_profile(input_file_path: str) -> rasterio.profiles.Profile:
    with rasterio.open(input_file_path) as raster:
        return raster.profile


def image_with_world_file_to_geotiff(input_file_path: str, output_file_path: str, input_crs: rasterio.crs.CRS = None,
                                     set_nodata_to_zero: bool = False):
    """
    Take an image with a world file and produce a GeoTiff with the same data and tagged with input_crs.
    This is useful when you have e.g. a jpg+jgw, and you know its CRS, because you can produce
    a single GeoTiff file with the proper CRS.
    If set_nodata_to_zero is True (False by default) the resulting GeoTiff will have the value 0 marked as nodata.
    """
    with rasterio.open(input_file_path) as raster:
        data = raster.read()
        # Same profile as the input file but with the new driver
        output_file_profile = raster.profile
        if set_nodata_to_zero:
            output_file_profile["nodata"] = 0
        output_file_profile["crs"] = input_crs
        output_file_profile["driver"] = "GTiff"
        with rasterio.open(output_file_path, 'w', **output_file_profile) as dst:
            dst.write(data)


def scale_raster_to_desired_res(input_file_path: str, output_file_path: str, desired_res_x: float,
                                desired_res_y: float = None):
    """
    Re-scale (resample) a raster file in input_file_path so it gets desired_res resolution (cell size).
    The resulting file will still be properly georeferenced, but the resolution of the cells will be different..
    If desired_res_y is omitted, desired_res_x will be used for both width and height
    """
    with rasterio.open(input_file_path) as raster:
        current_res_x = raster.profile["transform"][0]
        current_res_y = raster.profile["transform"][4]  # It will normally (always?) be a negative number
        scale_x = current_res_x / desired_res_x
        scale_y = -(current_res_y / desired_res_y)
        scale_raster(input_file_path, output_file_path, scale_x, scale_y)

def scale_raster(input_file_path: str, output_file_path: str, scale_factor_x: float, scale_factor_y: float = None,
                 resampling: rasterio.enums.Resampling = rasterio.enums.Resampling.nearest):
    """
    Re-scale (resample) a raster file in input_file_path by scale_factor_x and scale_factor_y and put the result
    in output_file_path. If scale_factor_y is not provided, scale_factor_x will be used for both width and height.
    If scale_factor < 0, the resulting file will have a lower resolution and vice versa.
    The resulting file will still be properly georeferenced, but the resolution of the cells will be different, and
    not necessarily the same in x and y, even if the original file had square cells.
    The resampling strategy is given by the parameter resampling (defaults to nearest neighbor)
    """
    if scale_factor_y is None:
        scale_factor_y = scale_factor_x

    with rasterio.open(input_file_path) as raster:

        logger.info(raster.profile)
        # Resample raster
        # Reading form a raster source into an array of a different size or with a specific out_shape,
        # is how we resample the data
        data = raster.read(
            out_shape=(
                raster.count,  # We keep the same number of bands of the input file
                # TODO: Make ceil/floor/even round? a parameter
                # Ceil: we do not lose any information, but for the latest row and column of pixels, we just
                # really had information about part of it in the original data so we are extrapolating ("inventing")
                # Floor: we may loose some information, but we don't have to extrapolate, all we do
                # is interpolate, so we do not "invent" information we didn't really have
                math.ceil(raster.height * scale_factor_y),
                math.ceil(raster.width * scale_factor_x)
            ),
            resampling=resampling
        )

        # data is a Numpy 2d ndarray (with more or less cells than the original after rescaling)
        logger.info(data.shape)

        # Transformation matrix. Without calculating this, the output file would be incorrectly georeferenced
        # because its profile would keep the original transform (that now is different)
        transform = raster.transform * raster.transform.scale(
            (raster.width / data.shape[-1]),
            (raster.height / data.shape[-2])
        )

        # Same profile as the input file but with the new transform
        output_file_profile = raster.profile
        output_file_profile["width"] = data.shape[-1]
        output_file_profile["height"] = data.shape[-2]
        output_file_profile["transform"] = transform
        logger.info(output_file_profile)

        with rasterio.open(output_file_path, 'w', **output_file_profile) as dst:
            dst.write(data)


def reproject_raster(input_file_path: str, output_file_path: str, dst_crs: rasterio.crs.CRS,
                     input_crs: rasterio.crs.CRS = None,
                     resampling: rasterio.enums.Resampling = rasterio.enums.Resampling.nearest,
                     src_nodata: Union[int, float] = None,
                     dst_nodata: Union[int, float] = None):

    # If input_crs is None it will read it from the input file (should be the most common option)
    with rasterio.open(input_file_path) as raster:
        if input_crs is None:
            input_crs = raster.profile["crs"]
        else:
            raster._crs = input_crs  # Yes, we are writing on a "private" attribute

        transform, width, height = rasterio.warp.calculate_default_transform(
            input_crs, dst_crs, raster.width, raster.height,
            *raster.bounds)

        if src_nodata is None:
            set_src_nodata = raster.nodata
        else:
            set_src_nodata = src_nodata

        if dst_nodata is None:
            set_dst_nodata = raster.nodata  # Index is band - 1 as it starts in 0.
        else:
            set_dst_nodata = dst_nodata

        # TODO: GeoTIFF just supports a single nodata value for all bands. If the input file has
        # a different nodata per band, the result will not be exactly as expected. Other output
        # format should be allowed (or at least a warning provided)
        kwargs = raster.meta.copy()
        kwargs.update({
            'driver': 'GTiff',  # Force writing a GeoTiff even if the original file is in another format
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height,
            'nodata': set_dst_nodata
        })

        with rasterio.open(output_file_path, 'w', **kwargs) as dst:
            for i in range(1, raster.count + 1):
                rasterio.warp.reproject(
                    source=rasterio.band(raster, i),
                    destination=rasterio.band(dst, i),
                    src_transform=raster.transform,
                    dst_transform=transform,
                    src_nodata=set_src_nodata,
                    dst_nodata=set_dst_nodata,
                    dst_crs=dst_crs,
                    resampling=resampling)


def get_geodesic_size_from_raster_profile(profile: rasterio.profiles.Profile, ellps:str = "WGS84") -> Tuple[float, float]:
    """
    Provides a rough approximation of the size in meters of a raster file with the given profile.
    It calculates the geodesic length of the line between the left-top and the right-bottom corners
    of the bbox of the given profile, and then uses its width and height to calculate an
    "average" cell resolution in meters (a single value, square cells).
    This is very rough. It uses the diagonal to avoid problems with global datasets where the width in the poles
    would be very small etc. If you need a more precise solution, **project your dataset** properly.
    Returns the length of the diagonal in meters, and the resolution in meters.
    """
    input_crs = profile["crs"]
    assert not input_crs.is_projected, "Use this function only to make a rough estimate of the dimensions of non-" \
                                       "projected raster datasets"
    geod = pyproj.Geod(ellps=ellps)
    left, top, right, bottom, current_res_x, current_res_y = get_bbox_from_raster_profile(profile)

    geodesic_diagonal = geod.line_length([left, right], [top, bottom])
    res = math.sqrt(geodesic_diagonal**2 / (profile["width"]**2 + profile["height"]**2))

    return geodesic_diagonal, res


# TODO: This is pretty slow, especially for large-ish datasets. If both have the same resolution
# it should be possible to use numpy to facilitate and make them faster the comparisons.
# Maybe there are some other optimizations for that and other cases.
def calculate_rmse(control_file_path: str, test_file_path: str, every_rows: int = 1,
                   every_cols: int = 1, band: int = 1) -> Tuple[float, float]:
    """
    Calculates de RMSE and BIAS between control_file_path (true values) and test_file_path
    (values we want to check). The test dataset could be the result of a predictive model, the result of some
    manipulation of the original which should keep its information (resampling, reprojecting...) etc.

    If every_cols and/or every_rows is larger than one, it will not test every cell, but every every_cols
    and every_rows cells. This may be sufficient for some purposes, and much faster.
    """
    band_idx = band - 1  # bands start in 1 for gdal/rasterio, but in the np array they start at 0

    with rasterio.open(control_file_path) as control:
        control_data = control.read()  # array with shape (bands, rows, columns)
        control_crs = control.profile["crs"]
        band_dtype = control.dtypes[band_idx]

        with rasterio.open(test_file_path) as test:
            test_data = test.read()
            test_crs = test.profile["crs"]

            # Without the always_xy=True this function fails for some combinations of control_crs and test_crs
            # (e.g., when the control is lat/lon and the test is rhealpix)
            transformer = pyproj.Transformer.from_crs(control_crs, test_crs, always_xy=True)

            sum_of_squares = 0
            sum_bias = 0
            n = 0

            logger.info(f"Nodata values for control: {control.nodatavals} and test: {test.nodatavals}")

            for row in range(0, control.height, every_rows):
                for col in range(0, control.width, every_cols):
                    control_value = control_data[band_idx, row, col]
                    if control_value != control.nodatavals[band_idx]:
                        # Get row,col in spatial coordinates
                        x, y = control.xy(row, col)

                        # Transform x,y to the crs of the test dataset
                        xtest, ytest = transformer.transform(x, y)
                        # Get the array index for that point in the test dataset
                        try:
                            r, c = test.index(xtest, ytest)
                            try:
                                test_value = test_data[band_idx, r, c]
                                if test_value != test.nodatavals[band_idx]:
                                    # If test_value and control_value are short integers, its difference may
                                    # produce an overflow (and thus incorrect values).
                                    if band_dtype in {'uint8', 'uint16', 'int8', 'int16'}:
                                        diff = int(test_value) - int(control_value)
                                    else:
                                        diff = test_value - control_value
                                    squared_diff = diff ** 2
                                    sum_of_squares += squared_diff
                                    sum_bias += diff
                                    n += 1
                            except IndexError as err:
                                logger.warning(f"x {x}, y {y}, xtest {xtest}, ytest {ytest}, r {r}, c {c}")
                                logger.warning(f"IndexError {err}. If this is occasional, This is just a rounding "
                                               f"issue or a, common, situation where after reprojecting to rhealpix, "
                                               f"the resulting dataset has lots of 'blank space' without a correspondence "
                                               f"with the original one.")
                        except ValueError:
                            logger.warning(f"Can't calculate row and col index for xtest {xtest}, ytest {ytest}")

            mean_sum_of_squares = sum_of_squares / n
            rmse = math.sqrt(mean_sum_of_squares)
            bias = sum_bias / n
            return rmse, bias