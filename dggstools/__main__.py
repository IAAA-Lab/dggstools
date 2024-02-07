from typing import Annotated, Optional

import rasterio
import typer

from rhealpixdggs.dggs import RHEALPixDGGS
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID

from dggstools.rhpx.vector_to_rhpx import vector_to_rhealpix, calculate_vector_raster_area_error
from dggstools.rhpx.raster_to_rhpx import raster_to_rhealpix, RescalingStrategy
from dggstools.rhpx.rhpxdataframes import RHEALPixDataFrameHelper
from dggstools.rhpx.utils.storage import geodataframe_to_geopackage, get_gpkg_rhpx_metadata

app = typer.Typer()


def _parse_rdggs(rdggs_str: str) -> RHEALPixDGGS:
    """
    Parse a string with the format:
    N_side/north_square/south_square
    Returns a RHEALPixDGGS object with the WG84 ellipsoid and the given parameters
    """
    n_side, north_square, south_square = rdggs_str.split("/")
    assert int(n_side) == 2 or int(n_side) == 3
    assert int(north_square) >= 0 and int(north_square) <= 3
    assert int(south_square) >= 0 and int(south_square) <= 3

    return RHEALPixDGGS(ellipsoid=WGS84_ELLIPSOID, north_square=int(north_square), south_square=int(south_square), N_side=int(n_side))


@app.command()
def vector_to_rhpx_raster(input_file_path: Annotated[str, typer.Argument()],
                          output_file_path: Annotated[str, typer.Argument()],
                          dst_resolution_idx: Annotated[int, typer.Argument()],
                          property_for_class: Annotated[Optional[str], typer.Argument()] = None,
                          fixed_value: Annotated[int, typer.Option()] = 1,
                          input_crs: Annotated[Optional[str], typer.Option()] = None,
                          layer: Annotated[Optional[str], typer.Option()] = None,
                          all_touched: Annotated[bool, typer.Option()] = False,
                          rdggs: Annotated[str, typer.Option()] = "3/1/0"):
    # layer is Union[str, int] but typer does not support Union
    # I will have to parse it manually

    # input_crs must be something that can be parsed by rasterio.crs.CRS.from_string or None
    if input_crs is not None:
        input_crs = rasterio.crs.CRS().from_string(input_crs)
    try:
        vector_to_rhealpix(_parse_rdggs(rdggs), input_file_path, output_file_path, dst_resolution_idx,
                           property_for_class, fixed_value, input_crs, layer, all_touched)
        result = "OK"
    except Exception as e:
        result = str(e)
    print(result)

@app.command()
def vector_raster_area_error(vector_file_path: Annotated[str, typer.Argument()],
                             raster_file_path: Annotated[str, typer.Argument()],
                             property_for_class: Annotated[Optional[str], typer.Argument()] = None,
                             fixed_value: Annotated[int, typer.Option()] = 1,
                             input_crs: Annotated[Optional[str], typer.Option()] = None,
                             band: Annotated[int, typer.Option()] = 1,
                             layer: Annotated[Optional[str], typer.Option()] = None):
    # layer is Union[str, int] but typer does not support Union
    # I will have to parse it manually

    # input_crs must be something that can be parsed by rasterio.crs.CRS.from_string or None
    if input_crs is not None:
        input_crs = rasterio.crs.CRS().from_string(input_crs)
    try:
        rmse, bias, total_vector_area, total_cells_area = calculate_vector_raster_area_error(vector_file_path, raster_file_path,
                                                                                             property_for_class, fixed_value,
                                                                                             input_crs,
                                                                                             band, layer)
        result = f"RMSE: {rmse} ({rmse/1000000} km²), bias: {bias} ({bias/1000000} km²), total_vector_area: {total_vector_area}, total_cells_area: {total_cells_area}"
    except Exception as e:
        result = str(e)
    print(result)


@app.command()
def raster_to_rhpx_raster(input_file_path: Annotated[str, typer.Argument()],
                          output_file_path: Annotated[str, typer.Argument()],
                          dst_resolution_idx: Annotated[int, typer.Argument()],
                          rescaling_strategy: Annotated[
                              RescalingStrategy, typer.Option()] = RescalingStrategy.TO_CLOSEST.value,
                          input_crs: Annotated[Optional[str], typer.Option()] = None,
                          resampling: Annotated[
                              rasterio.enums.Resampling, typer.Option()] = rasterio.enums.Resampling.nearest.value,
                          src_nodata: Annotated[Optional[int], typer.Option()] = None,
                          dst_nodata: Annotated[Optional[int], typer.Option()] = None,
                          rdggs: Annotated[str, typer.Option()] = "3/1/0"):

    # src_nodata and dst_nodata should be int | float, but typer does not support Union

    # input_crs must be something that can be parsed by rasterio.crs.CRS.from_string or None
    if input_crs is not None:
        input_crs = rasterio.crs.CRS().from_string(input_crs)

    try:
        raster_to_rhealpix(_parse_rdggs(rdggs), input_file_path, output_file_path, dst_resolution_idx,
                           rescaling_strategy, input_crs, resampling, src_nodata, dst_nodata)
        result = "OK"
    except Exception as e:
        result = str(e)
    print(result)

@app.command()
def raster_rhpx_to_geopackage(input_file_path: Annotated[str, typer.Argument()],
                              output_file_path: Annotated[str, typer.Argument()],
                              geo_id_column_name: Annotated[str, typer.Option()] = "cellid",
                              layer_name: Annotated[str, typer.Option()] = "data",
                              add_uid: Annotated[bool, typer.Option()] = False,
                              values_in_json: Annotated[bool, typer.Option()] = False,
                              store_nodata: Annotated[bool, typer.Option()] = False,
                              rdggs: Annotated[str, typer.Option()] = "3/1/0"):
    try:
        gdf = RHEALPixDataFrameHelper(_parse_rdggs(rdggs)).rhealpix_file_to_geodataframe(input_file_path,
                                                                                         geo_id_column_name, add_uid,
                                                                                         values_in_json, store_nodata)
        geodataframe_to_geopackage(gdf, output_file_path, layer_name)
        result = "OK"
    except Exception as e:
        result = str(e)
    print(result)

@app.command()
def print_gpkg_rhpx_metadata(input_file_path: Annotated[str, typer.Argument()]):
    try:
        metadata = get_gpkg_rhpx_metadata(input_file_path)
        result = str(metadata) + "\nOK"
    except Exception as e:
        result = str(e)
    print(result)



# TODO: Other possible commands

# def vector_to_rhpx_vector() # THIS IS JUSTA A BASIC VECTOR REPROJECTION, BUT A WRAPPER HERE WOULD SEEM OK
# def get_raster_rhpx_metadata(): THIS IS NOT YET IMPLEMENTED. NOT TRIVIAL, BUT NOT DIFFICULT. IT SHOULD PRODUCE THE SAME
# METADATA AS THE GEOPACKAGE COMMAND
# def vector_raster_line_error(): FAR LESS INTERESTING (AND TESTED) THAN THE AREA ERROR
# def rhealpix_grid_as_geodataframe(): # THERE ARE OTHER TOOLS TO GENERATE DGGS GRIDS. BUT THE FUNCTIONALITY IS AVAILABLE
# something for the auids


if __name__ == "__main__":
    app()

