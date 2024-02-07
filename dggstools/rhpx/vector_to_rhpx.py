from typing import Dict, Sequence

import fiona
import fiona.crs
import rasterio.crs
import shapely.geometry
from shapely.geometry import shape

from .raster_to_rhpx import *
from .rhpxutils import *
from .utils.vectorutils import *
from .utils.utils import *
from rasterio import features
from rasterio.warp import *

# TODO: ALLOW TO HAVE MORE THAN ONE PROPERTY AS VALUES FOR THE OUTPUT DATA (POSSIBLY AS DIFFERENT BANDS?)
# TODO: WHAT TO DO FOR POINTS AND LINES?

logger = logging.getLogger(__name__)


def __write_image(output_file_path: str, image: np.ndarray, input_crs: rasterio.crs.CRS, width: int, height: int,
                  transform: Affine):
    with rasterio.open(
            output_file_path, 'w',
            driver='GTiff',
            compress='DEFLATE',
            dtype=rasterio.uint16,  # For this function, uint8 could be enough, but this seems safer (at least if yo
            # don't know the number of feature classes in the vector file in advance
            crs=input_crs,
            transform=transform,
            count=1,
            width=width,
            height=height) as dst:
        dst.write(image, indexes=1)


def _get_feature_class(feature: Dict, property_for_value: str, fixed_value: int, props_to_integers: Dict[str, int])\
        -> Tuple[int, Dict[str, int]]:
    """
    Takes a Fiona feature and returns the value of a certain property:
    - The one named property_for_class, if that is not None and that property type is int.
    - fixed_value if property_for_class is None.
    - An integer stored in the dictionary props_to_integers for the value of property_for_class if it is a string.
      If that string is not in the dictionary, it adds it associated to a new int.

    Returns the value and the props_to_integers dictionary
    """
    if property_for_value is None:
        value = fixed_value
    else:
        try:
            value = int(feature["properties"][property_for_value])
        except ValueError:
            try:
                float(feature["properties"][property_for_value])
                raise TypeError(f"{property_for_value} is a float, but must be an int or str")
            except ValueError:
                # str is
                prop = feature["properties"][property_for_value]
                if prop not in props_to_integers:
                    props_to_integers[prop] = len(
                        props_to_integers.keys()) + 1  # Essentially this is an autoinc starting at 1
                value = props_to_integers[prop]
    return value, props_to_integers


def _geoms_and_values(input_features: List[Dict], property_for_value: str, fixed_value: int) -> \
        Tuple[List[Tuple[Dict, int]], Dict[str, int]]:
    """
    Takes input_features and returns a list with tuples with each geometry and
    - The value of property_for_class for that geometry, if it is an int
    - An int for each different value that property_for_class can take if it is a string
    - fixed_value if property_for_class is None
    """
    data = []
    props_to_integers = dict()  # Used only for non integer values
    for feature in input_features:
        value, props_to_integers = _get_feature_class(feature, property_for_value, fixed_value, props_to_integers)
        data.append((feature["geometry"], value))
    return data, props_to_integers


def vector_to_rhealpix(rdggs: RHEALPixDGGS, input_file_path: str, output_file_path: str, dst_resolution_idx: int,
                       property_for_class: str = None, fixed_value: int = 1, input_crs: rasterio.crs.CRS = None,
                       layer: Union[str, int] = None, all_touched: bool = False)\
        -> Tuple[int, Dict[str, int]]:
    """
    Take a vector file in input_file_path and produce a rasterized version in rhealpix in output_file_path:
    - With the RHEALPixDGGS defined in rdggs
    - With the rhealpix resolution number dst_resolution_idx
    - Using the field property_for_class to get the value used in the raster file. It will use fixed_value (default 1)
      if no property_for_class is provided.
    - Tries to get the input crs from the source data unless input_crs is provided.
    - With the layer fact_type or index given by layer in multilayer formats (e.g., GeoPackage). If none provide, the
    first one is used.
    - Uses an intermediate raster file which path will be derived from output_file_path adding _tmp before
      the extension.
    - The raster file produced uses uint16 as data type. This means that over 65535 classes it will not work as
      expected. This should be large enough for any sensible use of this function. uint8 would have been enough
      in many cases, in the future this could be parameterized for those situations when you know in advance the
      number of classes.
    - all_touched, by default False, means that only pixels whose center is within the polygon or that are selected
      by Bresenham’s line algorithm will be burned. That would be more accurate, but it is not always what you need.

    """
    rdggs_helper = RHEALPixDGGSHelper(rdggs)
    dst_resolution = rdggs.cell_width(dst_resolution_idx)

    with fiona.open(input_file_path, "r", layer=layer) as vectorfile:
        logger.info(vectorfile.bounds)
        input_features = [feature for feature in vectorfile]

        if input_crs is None:
            input_crs = rasterio.crs.CRS.from_string(fiona.crs.CRS.to_string(vectorfile.crs))
        logger.info("Input_crs: " + str(input_crs))

        left, top, right, bottom = bounds_to_left_top_right_bottom(vectorfile.bounds, input_crs)
        logger.info([left, top, right, bottom])

    # We calculate width and height on rhealpix when input data is not projected (or not in meters)
    # This way, we can calculate a resolution for the intermediate raster file which is close to the dst_resolution
    # (dst_resolution is in meters)
    if not input_crs.is_projected or not input_crs.linear_units == "metre":
        left_rheal, top_rheal, right_rheal, bottom_rheal = rdggs_helper.get_bbox_in_rhealpix(input_crs, left, top, right, bottom)
        width = round(abs(right_rheal - left_rheal) / dst_resolution)
        height = round(abs(top_rheal - bottom_rheal) / dst_resolution)
    else:
        width = round(abs(right - left) / dst_resolution)
        height = round(abs(top - bottom) / dst_resolution)


    # Just to prevent a division by zero, but if the calculated width or height are 0, something else is
    # probably wrong
    width = max(1, width)
    height = max(1, height)

    logger.info(f"w/h: {[width, height]}")
    logger.info(f"dst_resolution: {dst_resolution}")

    transform, width, height = calculate_default_transform(
       input_crs, input_crs, width, height,
       left=left, right=right, top=top, bottom=bottom)

    # Take geometries associated with a value derived from property_for_class or fixed_value
    data_to_include, props_to_integers = _geoms_and_values(input_features, property_for_class, fixed_value)

    image = rasterio.features.rasterize(
        data_to_include,
        transform=transform,
        all_touched=all_touched,
        out_shape=(height, width))

    __write_image(insert_suffix(output_file_path, "_tmp"), image, input_crs, width, height, transform)

    # TODO: area_or_point should be Pixel if the vector dataset is a point or line dataset (maybe?)
    # and Area if polygons (always?)

    # If input data is not projected, the rescaling strategy will fail, so we pass along dst_resolution_idx
    if not input_crs.is_projected or not input_crs.linear_units == "metre":
        res_idx = raster_to_rhealpix(rdggs, insert_suffix(output_file_path, "_tmp"), output_file_path,
                                     dst_resolution_idx=dst_resolution_idx)
    else:
        res_idx = raster_to_rhealpix(rdggs, insert_suffix(output_file_path, "_tmp"), output_file_path,
                                     rescaling_strategy=RescalingStrategy.TO_CLOSEST)
    return res_idx, props_to_integers


def vector_to_rhealpix_without_intermediate_raster(rdggs: RHEALPixDGGS, input_file_path: str, output_file_path: str, dst_resolution_idx: int,
                                                   property_for_class: str = None, fixed_value: int = 1, input_crs: rasterio.crs.CRS = None,
                                                   layer: Union[str, int] = None,
                                                   alignment_rhealpix_aware: bool = True, all_touched: bool = False) \
        -> Tuple[int, Dict[str, int]]:
    """
    DO NOT USE THIS FUNCTION WITHOUT READING CAREFULLY THIS.

    Rasterizes a vector file to rhealpix, but reprojecting the geometries and then burning those instead of
    via a intermediate raster file.
    In general, it is better to use the other function (vector_to_rhealpix):
    - The resulting datasets are different, but it is not clear which one is better. The differences are small and
      do not seem significant, at least in my tests so far.
    - This one is more fragile. Reprojecting geometries to rhealpix requires first that the original ones are
      topologically fit (they not always are) and then there are corner case (geometries that fall between
      two faces of the cube for example) that should be treated with caution.
    - This one is slower.
    - This one is far less tested and for now it should be considered unfinished work.
    - Do not mess with the alignment_rhealpix_aware default unless you are testing things. If you give it
      False value, the alignment will most certainly be wrong, at least in the vertical axis.
    - This one does not consider (yet) the area_or_point metadata of the output file
    - all_touched, by default False, means that only pixels whose center is within the polygon or that are selected
      by Bresenham’s line algorithm will be burned. That would be more accurate, but it is not always whay you need.

    - There is a possible exception: it works better for lineal features. If that is interesting for you, testing
    and improving further this function could be a good idea.
    """
    rdggs_helper = RHEALPixDGGSHelper(rdggs)
    dst_resolution = rdggs.cell_width(dst_resolution_idx)

    with fiona.open(input_file_path, "r", layer=layer) as vectorfile:
        logger.info(vectorfile.bounds)
        input_features = [feature for feature in vectorfile]

        if input_crs is None:
            input_crs = rasterio.crs.CRS.from_string(fiona.crs.CRS.to_string(vectorfile.crs))
        logger.info("Input_crs: " + str(input_crs))

        left, top, right, bottom = bounds_to_left_top_right_bottom(vectorfile.bounds, input_crs)
        logger.info([left, top, right, bottom])

    left_rheal, top_rheal, right_rheal, bottom_rheal = rdggs_helper.get_bbox_in_rhealpix(input_crs, left, top, right,
                                                                           bottom)
    width = round(abs(right_rheal - left_rheal) / dst_resolution)
    height = round(abs(top_rheal - bottom_rheal) / dst_resolution)

    # Just to prevent a division by zero, but if the calculated width or height are 0, something else is
    # probably wrong
    width = max (1, width)
    height = max(1, height)

    logger.info(f"w/h: {[width, height]}")
    logger.info(dst_resolution)

    dst_crs = rdggs_helper.rhealpixdef_to_proj_string()

    # Take geometries associated with a value derived from property_for_class or fixed_value
    data_to_include, props_to_integers = _geoms_and_values(input_features, property_for_class, fixed_value)
    data_to_include = [(rdggs_helper.project_and_clip_to_rhealpix(geom, input_crs.to_string()), value)
                       for geom, value in data_to_include]

    # Calculate the transform with these parameters. The idea is that the resulting dataset will
    # not only have the desired resolution (a cell size which corresponds to a given resolution in
    # the DGGS), but will also be aligned with the DGGS grid at that resolution (i.e., each cell in
    # the output raster corresponds exactly with a cell in the DGGS)
    transform, width, height = rasterio.warp.calculate_default_transform(
        input_crs, dst_crs, width, height,
        left=left, right=right, top=top, bottom=bottom,
        resolution=dst_resolution)
    # This new transform is not exactly aligned with the rhealpix grid. We align it:
    if alignment_rhealpix_aware:
        transform = rdggs_helper.align_transform(transform, dst_resolution_idx)
    else:
        # This results in wrong vertical alignments. Kept just for testing (and because I am not
        # really sure if it could be fixed)
        transform, width, height = rasterio.warp.aligned_target(transform, width, height, dst_resolution)

    # By default, only pixels whose center is within the polygon or that are selected by Bresenham’s line algorithm
    # will be burned in. It seems the most accurate option and it is the default (all_touched=False).
    image = rasterio.features.rasterize(
        data_to_include,
        transform=transform,
        all_touched=all_touched,
        out_shape=(height, width))

    __write_image(output_file_path, image, rasterio.crs.CRS.from_string(dst_crs), width, height, transform)

    return dst_resolution_idx, props_to_integers


# For now this might be used to initialize the optimal set of cuids for a given area in order to initialize an AUID
# for that area. However, is far slower than rasterizing that area at the highest desired resolution level,
# getting those CUIDS and then optimizing them before creating the AUID with them.
# I keep this because it might become the seed of a future variable resolution rasterization function which
# goes beyond simpliy rasterizing "plain" areas.
# Besides this, it is not very well designed/explained, and there is much room for improvement regarding performance,
# even while keeping the same algorithm. It might also have some errors yet.
def _vector_to_optimal_set_of_cuids(rdggs: RHEALPixDGGS, input_file_path: str, output_base_file_path: str,
                                    dst_resolution_idx: int,
                                    property_for_class: str,
                                    included_classes: Sequence[Union[int, str, float]],
                                    input_crs: rasterio.crs.CRS = None, layer: Union[str, int] = None) -> Sequence[str]:

    for i in range(0, dst_resolution_idx + 1):
        # A cell in resolution X can be NODATA and however be completely covered in
        # the highest resolution if we use all_touched = False, becase that rasterization
        # is more accurate (a cell is in or out if the border line includes its center).
        # We need all_touched = True for all the resolutions lower than the higher
        # resolution; this takes cells by excess, and this is what want.
        # If we don't do this, we may end up with less than optimal rasterziations, taht
        # do not include e.g. N1 but do include N11, N12, N13 and N14.
        # This happens, is not just a hypothetical case.
        if i == dst_resolution_idx:
            all_touched = False
        else:
            all_touched = True
        current_output_path = insert_suffix(output_base_file_path, f"-{i}")
        # This is relatively fast. We could gain something but reading the vector file just once
        # but most of the time is not spent here
        _, _ = vector_to_rhealpix(rdggs, input_file_path, current_output_path,
                                  input_crs=input_crs,
                                  dst_resolution_idx=i,
                                  property_for_class=property_for_class,
                                  layer=layer,
                                  all_touched=all_touched)

    highest_res_path = insert_suffix(output_base_file_path, f"-{dst_resolution_idx}")
    optimal_cuids = set()  # It won't have repetitions, so it should be faster than a list
    with rasterio.open(highest_res_path) as highest_resolution_raster:
        hr_left, hr_top, hr_right, hr_bottom, hr_res_x, hr_res_y = \
            get_bbox_from_raster_profile(highest_resolution_raster.profile)
        highest_resolution_raster_data = highest_resolution_raster.read()

        for current_resolution_idx in range(0, dst_resolution_idx + 1):
            current_res_path = insert_suffix(output_base_file_path, f"-{current_resolution_idx}")
            current_cell_side_in_highest_resolution_cells = rdggs.N_side ** (
                        dst_resolution_idx - current_resolution_idx)

            with rasterio.open(current_res_path) as current_resolution_raster:
                current_resolution_raster_data = current_resolution_raster.read()

                left, top, right, bottom, res_x, res_y = get_bbox_from_raster_profile(
                    current_resolution_raster.profile)

                # Different resolutions for the same vector file do not have the same extents. This is
                # so because we can't have "half DGGS cells" in the raster files. We need some padding
                # to make it easier to calculate which rows/cols in the highest resolution raster correspond
                # to ones in the lower resolution ones
                if left < hr_left:
                    hr_left_padding = round(abs(hr_left - left) / abs(hr_res_x))
                else:
                    hr_left_padding = -round(abs(left - hr_left) / abs(hr_res_x))
                if top < hr_top:
                    hr_top_padding = -round(abs(hr_top - top) / abs(hr_res_y))
                else:
                    hr_top_padding = round(abs(top - hr_top) / abs(hr_res_y))


                for row in range(0, current_resolution_raster.height):
                    for col in range(0, current_resolution_raster.width):

                        min_highest_row = row * current_cell_side_in_highest_resolution_cells
                        min_highest_col = col * current_cell_side_in_highest_resolution_cells

                        # We assume input raster has only one band (and hence the 0)
                        value = current_resolution_raster_data[0, row, col]
                        # logger.info(f"Value: {value}")

                        # if every value in the highest_resolution_raster is equal to value then
                        # we can keep the cell in the current_resolution_raster. But only if the number of
                        # highest resolution cells is equal to the maximum possible value, and all the
                        # values are equal among themselves
                        max_corresponding_cells = current_cell_side_in_highest_resolution_cells ** 2

                        def _calculate_corresponding_cells() -> int:
                            corresponding_cells = 0

                            for highest_row in range(min_highest_row - hr_top_padding,
                                                     min_highest_row - hr_top_padding +
                                                     current_cell_side_in_highest_resolution_cells):
                                for highest_col in range(min_highest_col - hr_left_padding,
                                                         min_highest_col - hr_left_padding +
                                                         current_cell_side_in_highest_resolution_cells):
                                    try:

                                        highest_value = highest_resolution_raster_data[
                                            0, highest_row, highest_col]
                                        if highest_value == value and highest_value in included_classes:
                                            # We don't really need to add up. We could always return
                                            # max_corresponding_cells
                                            # unless we return -1 (see the else). But this is cheap and for
                                            # logging/debugging has been useful
                                            corresponding_cells += 1
                                        else:
                                            return -1
                                    except IndexError:
                                        # It's fine, we are just in the padding
                                        pass
                            # logger.info(f"corresponding_cells: {corresponding_cells}, max_corresponding_cells: {max_corresponding_cells}")
                            return corresponding_cells

                        corresponding_cells = _calculate_corresponding_cells()
                        if corresponding_cells == max_corresponding_cells:
                            # TODO: This is quite inefficient. We are not just using the fact that the
                            # input data is aligned to the DGGS, so finding out which cellid corresponds to
                            # the current row and col could be done faster than this
                            x, y = current_resolution_raster.xy(row, col)
                            curr_cell = rdggs.cell_from_point(current_resolution_idx, (x, y))

                            # Inefficient, for now it's OK. Something based in regex could be faster.
                            add_curr_cell = True
                            for i in range(1, len(str(curr_cell))):
                                if str(curr_cell)[:-i] in optimal_cuids:
                                    add_curr_cell = False
                                    break

                            if add_curr_cell:
                                optimal_cuids.add(str(curr_cell))

    return sorted(optimal_cuids)


def calculate_vector_raster_area_error(vector_file_path: str, raster_file_path: str,
                                       property_for_class: str = None,
                                       fixed_value: int = 1,
                                       input_crs: rasterio.crs.CRS = None,
                                       band: int = 1,
                                       layer: Union[str, int] = None) -> Tuple[float, float, float, float]:
    """
    Take a vector file and a rasterized rhealpix file version and measures the area of
    each geometry in vector file and compares it with the areas of the cells which correspond to
    that geometry in the vector file.
    The correspondence is based on the value of property_for_class an layer in the vector file, and the
    given band (by default 1) in the raster file. This values identify feature classes, so the
    comparison essentially consists in measuring the difference in the areas of each feature
    class and aggregating all of them as rmse and as bias (assuming the vector file has the true value).
    Besides that, the total area of all the geometries in the vector file and all the cells with
    data in the raster file is returned, to provide a rough idea of how big are the rmse and bias.

    E.g. if a vector file has two feature classes, with values 1 and 2, this function calculates
    the square error between the total area of the geometries of class 1 (say 20 m²) and the
    total area of the  cells with value 1 in band 1 (say 21 m²), the same with class 2 (say 12 m² and
    15 m²), and then calculates  the rmse and bias properly adding those values.
    rmse: sqrt( ((21-20)² + (15-12)²) / 2 )
    bias: ((21-20) + (15-12)) / 2
    """
    band_idx = band - 1  # bands start in 1 for gdal/rasterio, but in the np array they start at 0
    vector_area = {}
    feature_classes = set()

    with fiona.open(vector_file_path, "r", layer=layer) as vectorfile:
        if input_crs is None:
            input_crs = rasterio.crs.CRS.from_string(fiona.crs.CRS.to_string(vectorfile.crs))

        assert (input_crs.is_projected and input_crs.linear_units == "metre") or not input_crs.is_projected, \
            f"If the crs of {vector_file_path} is projected, it must be in meters."

        props_to_integers = dict()  # Used only for non integer values
        for f in vectorfile:
            g = f["geometry"]
            assert g["type"] == "Polygon" or g["type"] == "MultiPolygon", \
            "Only polygons and multipolygons are allowed in this functions"

            feature_class, props_to_integers = _get_feature_class(f, property_for_class, fixed_value, props_to_integers)

            if feature_class not in feature_classes:
                feature_classes.add(feature_class)
                vector_area[feature_class] = 0
            if input_crs.is_projected:
                # If the projection is not an area preserving projection, this value could be
                # quite wrong: shape(g).area
                # So it is better that we take lon/lat coordinates and calculate its geodesic area
                with fiona.Env(OGR_ENABLE_PARTIAL_REPROJECTION="NO"):  # By default it is enabled
                    unprojected_g = fiona.transform.transform_geom(input_crs.to_string(),
                                                                   pyproj.CRS.from_epsg(4326).to_string(), g)
                if unprojected_g is not None:
                    # The next function defaults to WGS84, so I used 4326 in the previous statement
                    area, _ = get_geodesic_area_perimeter_of_geom(shape(unprojected_g))
                    vector_area[feature_class] += area
                else:
                    logger.info("Unprojection failed. Using the projected geometry to calculate the area.")
                    # If unprojecting has failed, we will use the projected area
                    vector_area[feature_class] += shape(g).area
            else:
                area, _ = get_geodesic_area_perimeter_of_geom(shape(g))
                vector_area[feature_class] += area

    with rasterio.open(raster_file_path) as raster:
        raster_crs = raster.profile["crs"]
        assert (raster_crs.is_projected and raster_crs.linear_units == "metre") or not raster_crs.is_projected, \
            f"If the crs of {raster_file_path} is projected, it must be in meters."
        _, _, _, _, res_x, res_y = get_bbox_from_raster_profile(raster.profile)
        if raster_crs.is_projected:
            cell_area = abs(res_x * res_y)
            try:
                if raster_crs.data["proj"] == "rhealpix":
                    # Correct cell_area with the areal distortion of rhealpix
                    cell_area = cell_area / RHEALPIX_MEAN_AREAL_DISTORTION
            except KeyError:
                # proj key is not in the data dictionary; it won't be rhealpix, we don't have to do anything
                pass
            # TODO: Other projections may have similar distortions. rhealpix is considered first as this is a
            # dggs library, but this area should be corrected (or calculated in another way) for other
            # crs which are not area-preserving
        else:
            # I am assuming WGS84, which is used by default in get_geodesic_size_from_raster_profile
            _, res = get_geodesic_size_from_raster_profile(raster.profile)
            cell_area = res * res
            logger.info("The area of the cells is an approximation because the raster dataset is not projected")
        raster_data = raster.read()  # array with shape (bands, rows, columns)

        num_cells = {v: np.count_nonzero(raster_data[band_idx] == v) for v in feature_classes}
        logger.info(f"Num cells with each value: {num_cells}")

        area_of_cells_with_values = {v: (np.count_nonzero(raster_data[band_idx] == v) * cell_area) for v in feature_classes}

    sum_of_squares = 0
    sum_bias = 0
    total_vector_area = 0
    total_cells_area = 0
    for fc in feature_classes:
        total_vector_area += vector_area[fc]
        total_cells_area += area_of_cells_with_values[fc]
        diff = area_of_cells_with_values[fc] - vector_area[fc]
        squared_diff = diff ** 2
        sum_of_squares += squared_diff
        sum_bias += diff

    mean_sum_of_squares = sum_of_squares / len(feature_classes)
    rmse = math.sqrt(mean_sum_of_squares)
    bias = sum_bias / len(feature_classes)
    return rmse, bias, total_vector_area, total_cells_area


def calculate_vector_raster_line_error(vector_file_path: str, raster_file_path: str,
                                       property_for_class: str = None,
                                       fixed_value: int = 1,
                                       input_crs: rasterio.crs.CRS = None,
                                       band: int = 1,
                                       every_feature: int = 1,
                                       layer: Union[str, int] = None) -> Tuple[float, float]:
    """
        Take a vector file with lineal features and a rasterized rhealpix file version and measures
        the average distance between each node in the features and the center of the raster cells
        where those nodes are.

        This function is pretty slow for large-ish datasets. If every_feature > 1, then only 1 in
        every_feature features will be taken into account. This provides an approximate error_msg and of
        course is much faster if every_feature is large enough.

        The correspondence is based on the value of property_for_class and layer in the vector file, and
        the given band (by default 1) in the raster file.

        If a node falls into an empty cell, or in a cell with a different feature class (probably because in the
        original vector file there where several features with different classes which were close to each
        other) that cell for now is not counted at all, but TODO it probably should add some error_msg.

        The result is the total sum of the distances between each node and its corresponding raster cell
        center divided by the num of nodes, and by the num of features.

        TODO: Cells in the raster that are marked as feature but where there is no feature are not counted
        as errors. I expect that a proper rasterization will minimize these. And it is not trivial to check
        for this, as looking only at the nodes of the lines would not be enough, and we would have to
        essentially rasterize the vector lines again... :-/
        """
    band_idx = band - 1  # bands start in 1 for gdal/rasterio, but in the np array they start at 0
    nodes = []  # to store every point in the input lines along with its feature class
    feature_classes = set()

    with fiona.open(vector_file_path, "r", layer=layer) as vectorfile:
        if input_crs is None:
            input_crs = rasterio.crs.CRS.from_string(fiona.crs.CRS.to_string(vectorfile.crs))

        assert (input_crs.is_projected and input_crs.linear_units == "metre") or not input_crs.is_projected, \
            f"If the crs of {vector_file_path} is projected, it must be in meters."

        props_to_integers = dict()  # Used only for non integer values
        skip_features_counter = every_feature
        total_nodes_tested = 0
        total_features_tested = 0
        for f in vectorfile:
            if skip_features_counter < every_feature:
                skip_features_counter += 1
            else:
                skip_features_counter = 1
                total_features_tested +=1
                g = f["geometry"]
                assert g["type"] == "LineString" or g["type"] == "MultiLineString", \
                "Only lines and multilines are allowed in this functions"

                feature_class, props_to_integers = _get_feature_class(f, property_for_class, fixed_value, props_to_integers)

                if feature_class not in feature_classes:
                    feature_classes.add(feature_class)
                if input_crs.is_projected:
                    # No projection preserves distances perfectly (though some do better than others).
                    # So it is safer that we take lon/lat coordinates
                    with fiona.Env(OGR_ENABLE_PARTIAL_REPROJECTION="NO"):  # By default it is enabled
                        g = fiona.transform.transform_geom(input_crs.to_string(),
                                                       pyproj.CRS.from_epsg(4326).to_string(), g)
                # Notice that I am assuming WGS84 in case it is not projected

                line_or_multiline = shapely.geometry.shape(g)
                for point in line_or_multiline.coords:
                    total_nodes_tested += 1
                    nodes.append((point, feature_class))

    with rasterio.open(raster_file_path) as raster:
        raster_crs = raster.profile["crs"]
        assert (raster_crs.is_projected and raster_crs.linear_units == "metre") or not raster_crs.is_projected, \
            f"If the crs of {raster_file_path} is projected, it must be in meters."
        _, _, _, _, res_x, res_y = get_bbox_from_raster_profile(raster.profile)

        raster_data = raster.read()  # array with shape (bands, rows, columns)
        if raster_crs.is_projected:
            cell_width = res_x
            transformer = pyproj.Transformer.from_crs(pyproj.CRS.from_epsg(4326), raster_crs, always_xy=True)
            inv_transformer = pyproj.Transformer.from_crs(raster_crs, pyproj.CRS.from_epsg(4326), always_xy=True)
        else:
            _, cell_width = get_geodesic_size_from_raster_profile(raster.profile)

        sum_of_distances_per_fc = {}
        total_distance = 0
        for p, fc in nodes:
            # p is in EPSG:4326 (we have transformed it before)
            # transform to the raster_crs if necessary (i.e., if raster is projected):
            if raster_crs.is_projected:
                px, py = transformer.transform(p[0], p[1])
            else:
                px = p[0]
                py = p[1]

            # row and col of the cell corresponding to the feature node point
            r, c = raster.index(px, py)  # TODO X AND Y COULD BE INVERTED?
            # Center of the r,c cell (it will normally not be p[0],p[1]. The difference is the error_msg we measure
            x, y = raster.xy(r, c)
            if raster_crs.is_projected:
                # We want x,y (center of the cell) in EPSG:4326, no in the raster crs
                x, y = inv_transformer.transform(x, y)

            try:
                raster_value = raster_data[band_idx, r, c]
                if raster_value == fc:
                    dist = get_geodesic_distance(x, y, p[0], p[1])
                else:
                    # The node falls in a raster cell which is not identified with the same
                    # feature class (probably it falls in a NODATA cell). This is also an error_msg,
                    # but we can't measure it with the distance to the center of that cell. This is frequent
                    # enough to deserve consideration.
                    # To account for it, I will add the a cell width as an error_msg
                    dist = cell_width
                total_distance += dist
                if fc not in sum_of_distances_per_fc:
                    sum_of_distances_per_fc[fc] = dist
                else:
                    sum_of_distances_per_fc[fc] += dist

            except IndexError as err:
                logger.warning(f"IndexError {err}. If this is occasional, I assume, for now, that it is a rounding issue.")

        mean_distance_per_node = total_distance / total_nodes_tested
        mean_distance_per_feature = total_distance / total_features_tested
        return mean_distance_per_node, mean_distance_per_feature
