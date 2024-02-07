import os
from enum import Enum

import rasterio.crs

from .rhpxutils import *
from .utils.rasterutils import *

logger = logging.getLogger(__name__)


# TODO: Generate two rasters (or something like that), one with the data reprojected to rhealpix and the other
# one with the error_msg per cell (information lost, extrapolated...) could be interesting. This is something
# that needs to be carefully thought

# TODO: Think about the implications of the different resampling strategies, as the differences for
# different datasets/domains are clear. This must be clearly communicated to the users of the
# reprojected/dggs data, specially if the data was downsampled (to lower resolution than the original).


class RescalingStrategy(Enum):
    TO_HIGHER = 1
    TO_LOWER = 2
    TO_CLOSEST = 3  # It will be the higher or the lower, the one which is closest


def _reproject_raster_to_rhealpix(rdggs: RHEALPixDGGS, input_file_path: str, output_file_path: str,
                                  dst_resolution_idx: int, input_crs: rasterio.crs.CRS = None,
                                  resampling: rasterio.enums.Resampling = rasterio.enums.Resampling.nearest,
                                  src_nodata: Union[int, float] = None, dst_nodata: Union[int, float] = None):
    """
    Reprojects a raster file in input_file_path to a raster file in output_file_path that will be:
    - Reprojected to rdggs (an rHEALPix based dggs). If input_crs is None, it will try to get the CRS from
    the input file (as long as it is defined there, this should be the most common option).
    - With the rHEALPix resolution which index is dst_resolution_idx (as determined by rdggs.cell_width())
    as the cell size.
    - In the GeoTiff format.
    """
    rdggs_helper = RHEALPixDGGSHelper(rdggs)
    dst_resolution = rdggs.cell_width(dst_resolution_idx)

    # If input_crs is None it will read it from the input file (should be the most common option)
    if input_crs is None:
        open_as = "r"
    else:
        # If input_crs is not None we must modify the input file to give it a CRS, so we **must**
        # open it in "r+" mode
        open_as = "r+"

    with rasterio.open(input_file_path, open_as) as raster:
        # Take the input file top, left, bottom and right coordinates of the input raster
        left, top, right, bottom, res_x, res_y = get_bbox_from_raster_profile(raster.profile)

        logger.info(f"BBOX OF INPUT FILE {[left, top, right, bottom, res_x, res_y]}")
        logger.info(f"PROFILE OF INPUT FILE {raster.profile}")

        if input_crs is None:
            input_crs = raster.profile["crs"]
        else:
            # The input file does not have a CRS or we want to change it. We set it or the rest of the
            # method will fail
            raster.crs = input_crs

        logger.info(f"INPUT CRS {input_crs}")
        dst_crs = rasterio.crs.CRS.from_string(rdggs_helper.rhealpixdef_to_proj_string())

        # A bit ad-hoc, and I don't have this problem with -180, 90 or -90 (or with some datasets with
        # global coverage and values much closer to 180 than the one I have to use here). For now this fixes
        # a problem with some datasets with some large extents
        if not input_crs.is_projected:
            max_right = 180.0 - 1e-1  # Depending on the dst_resolution, it seems, this could be a bit larger (not much)
            if right >= max_right:
                logger.warning(f"The right of the bbox from {input_file_path} is {right}; we change it to {max_right} "
                               f"to prevent a failure when reprojecting. You might want to double check the results.")
                right = max_right

        # Calculate the new transform
        transform, width, height = rasterio.warp.calculate_default_transform(
            input_crs, dst_crs, raster.width, raster.height,
            left=left, right=right, top=top, bottom=bottom,
            resolution=dst_resolution)

        # This new transform is not aligned with the rhealpix grid. We align it:
        transform = rdggs_helper.align_transform(transform, dst_resolution_idx)


        logger.info(f"After alignment: {[left, top, right, bottom]}")
        logger.info(f"After alignment: {transform}, {[width, height, dst_resolution]}")

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
            'compress': 'DEFLATE',
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
            # TOFIX: AREA_OR_POINT DOES NOT GET WRITTEN TO THE OUTPUT FILE. But I can write anything else
            # I haven't found a way to write this in dst.
            # TODO: Once I do, parameterize this function (and raster_to_rhealpix and maybe others) so the
            # client code can decide which one is used
            # None of this, among others, work
            # dst.update_tags(AREA_OR_POINT="Point")
            # dst.update_tags(**{'AREA_OR_POINT': "Point"})
            # dst.update_tags(RasterPixelIsPoint='<MDI "AREA_OR_POINT">Point</MDI>')
            # dst.update_tags(TIFFTAG_GDL_METADATA="Point")
            # dst.update_tags(GTRasterTypeGeoKey=1)
            # I suppose that writing that metadta with gdal would be the solution, but I don't want to
            # use gdal until I have a good solution for deployment (e.g. a dockerfile)
            # It could also be a bug in GDAL. The version I use (what comes with Ubuntu 20.04) is outdated
            # and there seem to be some bugfixes related to these issues




# TODO: If the original file is not a geotiff, we should verify that the output file has the tif extension or adding it
def raster_to_rhealpix(rdggs: RHEALPixDGGS, input_file_path: str, output_file_path: str, dst_resolution_idx: int = -1,
                       rescaling_strategy: RescalingStrategy = RescalingStrategy.TO_CLOSEST,
                       input_crs: rasterio.crs.CRS = None,
                       resampling: rasterio.enums.Resampling = rasterio.enums.Resampling.nearest,
                       src_nodata: Union[int, float] = None, dst_nodata: Union[int, float] = None) -> int:
    """
    Transform a raster file in input_file_path to a raster file in output_file_path:
    - Reprojected to rdggs (an rHEALPix based dggs). If input_crs is None, it will try to get the CRS from
    the input file (as long as it is defined there, this should be the most common option).
    - With the resolution fixed to one of the rHEALPix ones. By default the closest resolution to that
    of the original file, but you can choose a different rescaling_strategy, or even provide a specific
    resolution (dst_resolution_idx). If you provide a dst_resolution_idx, the rescaling_strategy
    argument will be ignored.
    - With the chosen resampling strategy (by default nearest neighbour).
    - In the Geotiff format.
    - Uses the src_nodata and dst_nodata values as the nodata for the input/output files; if any of them is None,
      it will try to read it from the input file and use that for the output. This is the best option: if your
      input file does not have the proper nodata value set in its metadata, you should try to fix it in your file
      and using these parameters just as last resort

     Returns the resolution used for the output dggs file index
    """
    if dst_resolution_idx == -1:  # If we don't get a specific resolution, we calculate one following rescaling_strategy
        rdggs_helper = RHEALPixDGGSHelper(rdggs)
        profile = get_raster_profile(input_file_path)
        logger.info(f"Input file profile: {profile}")
        left, top, right, bottom, current_res_x, current_res_y = get_bbox_from_raster_profile(profile)
        input_crs = profile["crs"]

        assert (input_crs.is_projected and input_crs.linear_units == "metre") or not input_crs.is_projected,\
            f"If the crs of {input_file_path} is projected, it must be in meters. Alternatively, you can provide " \
            f"an explicit value for dst_resolution_idx."

        if input_crs.is_projected and input_crs.linear_units == "metre":
            # For square cells it is very simple, but if current_res_x and abs(current_res_y) are different
            # we must choose one for the cells of rehalpix (which are squares). Choosing the best (smallest) one
            # is optimal from the point of view of keeping the information in the original file
            best_resolution = min(current_res_x, abs(current_res_y))
            # We could also take the average or follow other strategies...
            # avg_resolution = (current_res_x + abs(current_res_y)) / 2
        elif not input_crs.is_projected:
            # For non projected data, we can make a rough calculation of the resolution in meters, but it would
            # be better to project the data properly before using this function
            logger.info(f"Roughly estimating a resolution for this non-projected dataset. It would be better to "
                        f"project it")
            _, best_resolution = get_geodesic_size_from_raster_profile(profile)

        if rescaling_strategy == RescalingStrategy.TO_HIGHER:
            dst_resolution_idx, dst_res = rdggs_helper.get_closest_higher_resolution(best_resolution)
        elif rescaling_strategy == RescalingStrategy.TO_LOWER:
            dst_resolution_idx, dst_res = rdggs_helper.get_closest_lower_resolution(best_resolution)
        elif rescaling_strategy == RescalingStrategy.TO_CLOSEST:
            dst_resolution_idx, dst_res = rdggs_helper.get_closest_resolution(best_resolution)
        else:
            raise ValueError("Unknown rescaling strategy")

        logger.info(f"DGGS_RES={dst_res}")

    _reproject_raster_to_rhealpix(rdggs, input_file_path, output_file_path, dst_resolution_idx, input_crs, resampling,
                                  src_nodata, dst_nodata)
    return dst_resolution_idx

