import json
import os

import geopandas
import pyproj
import sqlalchemy
import sqlite3

from geopandas import GeoDataFrame
from typing import List, Dict, Any

import dggstools.rhpx.rhpxutils as rhpxutils
from dggstools.rhpx.rhpxdataframes import RHEALPixDataFrameHelper
from dggstools.rhpx.utils.rasterutils import get_raster_profile


def geodataframe_to_postgis(gdf: GeoDataFrame, table_name: str, username: str, password: str, host: str, port:int,
                            database: str, if_exists: str="fail", schema: str = "public", index:bool = False,
                            index_label: List[str] = None, chunksize: int = None,
                            dtype: Dict[str, sqlalchemy.types.TypeEngine] = None):
    """
       Saves gdf as table_name in a PostGIS DB.
       It stores the gdf "as is" (using GeoPandas to_postgis function).
       TODO: IT FAILS IF YOU TRY TO USE A GDF WITH A NON-EPSG CRS (E.G., RHEALPIX DEFINED AS A WKT STRING) UNLESS
       gdf.crs is previously set to None.
       """
    engine = sqlalchemy.create_engine(f"postgresql://{username}:{password}@{host}:{port}/{database}")
    gdf.to_postgis(table_name, engine, if_exists=if_exists, schema=schema, index=index, index_label=index_label,
                   chunksize=chunksize, dtype=dtype)


def geodataframe_to_geopackage(gdf: GeoDataFrame, output_file_path: str, layer_name:str ="data"):
    """
    Saves gdf as a GeoPackage in output_file_path, with the data in table named layer_name.
    Before that, it deletes output_file_path if it exists.
    Stores de gdf "as is" (using GeoPandas to_file function), but it also includes
    gdf.attrs in the geopackage_metadata table, following (more or less) the GeoPackage
    standard <http://www.geopackage.org/guidance/extensions/metadata.html>.
    """
    try:
        os.remove(output_file_path)
    except:
        # It should be OK. It does not exist, so we can't remove it before creating it again.
        # And if it is a different problem (e.g., it is a directory, you lack permissions...)
        # The problem will arise soon after this
        pass

    gdf.to_file(output_file_path, mode='w', driver='GPKG', layer=layer_name)

    # gdf.attrs don't get written automatically, so we have to add them
    con = sqlite3.connect(output_file_path)
    with con:
        con.execute(f"insert into gpkg_extensions (table_name, extension_name, definition, scope) "
                    f"values ('gpkg_metadata', 'gpkg_metadata', "
                    f"'http://www.geopackage.org/spec120/#extension_metadata', "
                    f"'read-write');")
        con.execute(f"insert into gpkg_extensions (table_name, extension_name, definition, scope) "
                    f"values ('gpkg_metadata_reference', 'gpkg_metadata', "
                    f"'http://www.geopackage.org/spec120/#extension_metadata', "
                    f"'read-write');")
        con.execute(f"CREATE TABLE gpkg_metadata (id INTEGER PRIMARY KEY AUTOINCREMENT,"
                    f"md_scope TEXT NOT NULL DEFAULT 'dataset', md_standard_uri TEXT NOT NULL, "
                    f"mime_type TEXT NOT NULL DEFAULT 'text/xml',  metadata TEXT NOT NULL DEFAULT '');")
        con.execute(f"insert into gpkg_metadata (id, md_scope, md_standard_uri, mime_type, metadata) "
                    f"values (1,'dataset','http://www.iaaa.es/geo2dggs/spec/1','application/json', ?);",
                    [json.dumps(gdf.attrs)])

    con.close()


def rhealpix_to_geopackage(input_file_path: str, output_file_path: str, geo_id_column_name: str = "cellid",
                           layer_name: str ="data",
                           add_uid: bool = False, values_in_json: bool = False, store_nodata: bool = False):
    profile = get_raster_profile(input_file_path)
    input_crs = profile["crs"]
    attrs = rhpxutils.get_gdf_attrs_from_rhealpix_file(input_file_path)
    input_dggs = rhpxutils.pyproj_crs_to_rdggs(pyproj.CRS(input_crs), int(attrs["rhealpixdggs"]["n_side"]))

    # Load the rhealpix raster file into a geodataframe
    gdf = RHEALPixDataFrameHelper(input_dggs).rhealpix_file_to_geodataframe(input_file_path, geo_id_column_name, add_uid,
                                                                            values_in_json, store_nodata)
    # # Save the geodataframe to a geopackage
    geodataframe_to_geopackage(gdf, output_file_path, layer_name)


def geopackage_to_rhealpix(input_file_path: str, output_file_path: str, nodata: int | float = 0):
    gdf = geopandas.read_file(input_file_path, engine="pyogrio")
    # # Add the attrs (metadata) to the geodataframe
    rhpx_metadata = get_gpkg_rhpx_metadata(input_file_path)
    gdf.attrs.update(rhpx_metadata)
    rdggs = rhpxutils.gdf_attrs_to_rdggs(rhpx_metadata)

    # Save the geodataframe to a rhealpix raster file
    RHEALPixDataFrameHelper(rdggs).geodataframe_to_rhealpix_file(gdf,
                                                                 output_file_path,
                                                                 metadata_dict=rhpx_metadata,
                                                                 nodata=nodata)


def get_gpkg_rhpx_metadata(input_file_path: str) -> dict[str, Any]:
    """
    Returns a dictionary with the rhealpix metadata (gdf attrs) of a geopackage that was saved with
    geodataframe_to_geopackage.
    """
    con = sqlite3.connect(input_file_path)
    with con:
        cur = con.cursor()
        cur.execute(f"select metadata from gpkg_metadata where id=1;")
        result = cur.fetchone()
        metadata_json = result[0]

    con.close()
    return json.loads(metadata_json)