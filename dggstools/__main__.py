from typing import Annotated, Optional

import rasterio
import typer

from rhealpixdggs.dggs import RHEALPixDGGS
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID

from dggstools.rhpx.rhpxutils import get_gdf_attrs_from_rhealpix_file
from dggstools.rhpx.vector_to_rhpx import vector_to_rhealpix, calculate_vector_raster_area_error
from dggstools.rhpx.raster_to_rhpx import raster_to_rhealpix, RescalingStrategy
from dggstools.rhpx.utils.storage import  get_gpkg_rhpx_metadata, rhealpix_to_geopackage, \
    geopackage_to_rhealpix

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
def vec_to_rhpx_ras(input_file_path: Annotated[str, typer.Argument()],
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
def vec_ras_area_error(vector_file_path: Annotated[str, typer.Argument()],
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
def ras_to_rhpx_ras(input_file_path: Annotated[str, typer.Argument()],
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
def ras_rhpx_to_vec_rhpx(input_file_path: Annotated[str, typer.Argument()],
                         output_file_path: Annotated[str, typer.Argument()],
                         geo_id_column_name: Annotated[str, typer.Option()] = "cellid",
                         layer_name: Annotated[str, typer.Option()] = "data",
                         add_uid: Annotated[bool, typer.Option()] = False,
                         values_in_json: Annotated[bool, typer.Option()] = False,
                         store_nodata: Annotated[bool, typer.Option()] = False):
    try:
        rhealpix_to_geopackage(input_file_path,
                               output_file_path, geo_id_column_name, layer_name, add_uid, values_in_json, store_nodata)
        result = "OK"
    except Exception as e:
        result = str(e)
    print(result)

@app.command()
def vec_rhpx_to_ras_rhpx(input_file_path: Annotated[str, typer.Argument()],
                         output_file_path: Annotated[str, typer.Argument()],
                         nodata: Annotated[float, typer.Option()] = 0):
    # nodata can be int or float, but typer does not support Union types
    try:
        geopackage_to_rhealpix(input_file_path, output_file_path, nodata)
        result = "OK"
    except Exception as e:
        result = str(e)
    print(result)


@app.command()
def print_vec_rhpx_metadata(input_file_path: Annotated[str, typer.Argument()]):
    try:
        metadata = get_gpkg_rhpx_metadata(input_file_path)
        result = str(metadata) + "\nOK"
    except Exception as e:
        result = str(e)
    print(result)

@app.command()
def print_ras_rhpx_metadata(input_file_path: Annotated[str, typer.Argument()]):
    try:
        metadata = get_gdf_attrs_from_rhealpix_file(input_file_path)
        result = str(metadata) + "\nOK"
    except Exception as e:
        result = str(e)
    print(result)


if __name__ == "__main__":
    app()

