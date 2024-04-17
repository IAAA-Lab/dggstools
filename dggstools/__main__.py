from enum import Enum
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

INPUT_CRS_HELP = "If not provided, the input file CRS will be used. If provided it must be a string parseable by rasterio.crs.CRS.from_string()."
RDGGS_HELP = "rHEALPix system parameters: n_side/north_square/south_square"

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
                    dst_resolution_idx: Annotated[int, typer.Argument(min=0)],
                    property_for_class: Annotated[Optional[str], typer.Option(help="Name of the attribute in the vector file that will be used for the values in the output raster file. If this parameter is not present, the --fixed-value will be used instead.")] = None,
                    fixed_value: Annotated[int, typer.Option(help="If --property-for-class is not used, then this is the value that will be used for the cells that correspond to any polygon in the input file.")] = 1,
                    input_crs: Annotated[Optional[str], typer.Option(help=INPUT_CRS_HELP)] = None,
                    layer: Annotated[Optional[str], typer.Option(help="If the input file is multilayer, e.g. a GeoPackage, name or index of the layer we want. By default, the first one is used.")] = None,
                    rdggs: Annotated[str, typer.Option(help=RDGGS_HELP)] = "3/1/0"):
    """
    Transforms a vector dataset with polygons in a common GIS format and reference system to a rHEALPix GeoTIFF.
    This GeoTIFF has a rasterization of the polygons following the constraints of the rHEALPix system (projection,
    valid resolution and grid alignment).

    --dst_resolution_idx is used to choose the resolution of the output file, by specifying
    the resolution index (0 is the lowest) we want.

    """
    # layer is Union[str, int] but typer does not support Union
    # I will have to parse it manually

    # input_crs must be something that can be parsed by rasterio.crs.CRS.from_string or None
    if input_crs is not None:
        input_crs = rasterio.crs.CRS().from_string(input_crs)

    if layer is not None:
        # typer does not support Union(str,int) as the type, so we have to do a manual parsing
        try:
            layer = int(layer)
        except ValueError:
            # It is OK. It is not an int, but a str is also a valid type
            pass

    try:
        vector_to_rhealpix(_parse_rdggs(rdggs), input_file_path, output_file_path, dst_resolution_idx,
                           property_for_class, fixed_value, input_crs, layer)
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
    """
    Takes a vector file and a rasterized rHEALPix version (as produced by the vec-to-rhpx-ras command) and:
     - measures the area of each geometry in vector file;
     - compares each of these areas with the areas of the cells which correspond to that geometry in the vector file.
    This is an experimental, not thoroughly tested and barely documented command, and
    should be used just for testing purposes.

    It will work for raster files which are not in rHEALPix too, but the results may be less precise depending on the
    areal distortion of their projection.
    """
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
                    rescaling_strategy: Annotated[
                              RescalingStrategy, typer.Argument()] = RescalingStrategy.TO_CLOSEST.value,
                    dst_resolution_idx: Annotated[Optional[int], typer.Option(min=0, help="If not provided, the rescaling-strategy will be used instead.")] = None,
                    input_crs: Annotated[Optional[str], typer.Option(help=INPUT_CRS_HELP)] = None,
                    # resampling is passed as a string, because rasterio.enums-Resampling is an IntEnum (not a StrEnum) and using numbers is less user-friendly
                    resampling: Annotated[
                        str, typer.Option(help=f"{[r.name for r in rasterio.enums.Resampling]}")] = rasterio.enums.Resampling.nearest.name,
                    rdggs: Annotated[str, typer.Option(help=RDGGS_HELP)] = "3/1/0"):
    """
    Transforms a raster dataset in a common GIS format and reference system to a rHEALPix GeoTIFF. This includes: warping to the rHEALPix projection, resampling to one of the allowed rHEALPix
    resolutions (which depend on the rHEALPix system being used) and aligning to that rHEALPix grid.

    --rescaling-strategy allows to choose which resolution to use in the output file among the ones
    allowed by the rHEALPix system: the default is to use the closest one (higher or lower), but we can
    choose the closest higher or the nex closest lower resolution.
    --dst_resolution_idx is an alternative way to choose the resolution of the output file, by specifying
    the resolution index (0 is the lowest) we want. This is optional and by default is a -1, which means
    that the --rescaling-strategy will be used instead.

    The resulting GeoTIFF includes metadata about the rHEALPix system it uses, but it is still
    a normal GeoTIFF: tools which can't read those metadata can use it as any other GeoTIFF,
    as long as they support the rHEALPix projection.
    """

    # src_nodata and dst_nodata should be int | float, but typer does not support Union

    # input_crs must be something that can be parsed by rasterio.crs.CRS.from_string or None
    if input_crs is not None:
        input_crs = rasterio.crs.CRS().from_string(input_crs)

    if dst_resolution_idx is None:
        dst_resolution_idx = -1

    try:
        raster_to_rhealpix(_parse_rdggs(rdggs), input_file_path, output_file_path, dst_resolution_idx,
                           rescaling_strategy, input_crs, rasterio.enums.Resampling[resampling])
        result = "OK"
    except Exception as e:
        result = str(e)
    print(result)

@app.command()
def ras_rhpx_to_vec_rhpx(input_file_path: Annotated[str, typer.Argument()],
                         output_file_path: Annotated[str, typer.Argument()],
                         geo_id_column_name: Annotated[str, typer.Option(help="Name of the column that will contain the cell identifier.")] = "cellid",
                         layer_name: Annotated[str, typer.Option(help="Name of the table in the output GeoPackage.")] = "data",
                         add_uid: Annotated[bool, typer.Option(help="If True, adds a column named uuid with a random UUID4 for each row in the output.")] = False,
                         values_in_json: Annotated[bool, typer.Option(help="If True, all bands are put together in a single column, in JSON format.")] = False):
    """
    Transforms a rHEALPix GeoTIFF dataset produced by dggstools to a vector dataset in the GeoPackage format.

    This GeoPackage has a point for each cell in the GeoTIFF (its centroid) and is equivalent to the raster
    version, up to the point that you can obtain back the original raster dataset from this GeoPackage, with
    only minor differences due to implementation details of the GeoTIFF format.

    Besides this, the GeoPackage is a normal vector dataset, that can be processed by any GIS application, as long as
    it can use the rHEALPix projection.

    The GeoPackage will have a column with the identifier (address) for each cell (--geo-id-column-name),
    and additional columns corresponding to the bands of the GeoTIFF (unless --values-in-json is used, in which case
    all bands will be in a single column named all_bands, in JSON format).
    """
    try:
        rhealpix_to_geopackage(input_file_path,
                               output_file_path, geo_id_column_name, layer_name, add_uid, values_in_json)
        result = "OK"
    except Exception as e:
        result = str(e)
    print(result)

@app.command()
def vec_rhpx_to_ras_rhpx(input_file_path: Annotated[str, typer.Argument()],
                         output_file_path: Annotated[str, typer.Argument()],
                         nodata: Annotated[float, typer.Option(help="Value for the NODATA pixels in the output raster; if you want it to be used as an INTEGER, don't write a decimal point.")] = 0.0):
    """
    Transforms a vector dataset in rHEALPix produced by dggstools with the ras-rhpx-to-vec-rhpx command, to
    a raster dataset in GeoTIFF which is very similar to the one that was used as the original input to that operation.
    """
    # nodata can be int or float, but typer does not support Union types, so:
    if isinstance(nodata, int):
        nodata = int(nodata)

    try:
        geopackage_to_rhealpix(input_file_path, output_file_path, nodata)
        result = "OK"
    except Exception as e:
        result = str(e)
    print(result)


@app.command()
def print_vec_rhpx_metadata(input_file_path: Annotated[str, typer.Argument()]):
    """
    Takes a GeoPackage file produced by dggstools and prints the metadata that dggstools stores in it. This metadata
    are necessary to store some rHEALPix system specific information, and some other information that can be useful
    if you want the original raster file back.
    """
    try:
        metadata = get_gpkg_rhpx_metadata(input_file_path)
        result = str(metadata) + "\nOK"
    except Exception as e:
        result = str(e)
    print(result)

@app.command()
def print_ras_rhpx_metadata(input_file_path: Annotated[str, typer.Argument()]):
    """
    Takes a GeoTIFF file produced by dggstools and prints the metadata that dggstools stores in it. This metadata
    are necessary to store some rHEALPix system specific information.
    """
    try:
        metadata = get_gdf_attrs_from_rhealpix_file(input_file_path)
        result = str(metadata) + "\nOK"
    except Exception as e:
        result = str(e)
    print(result)


if __name__ == "__main__":
    app()

