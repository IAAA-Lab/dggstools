import logging

import geopandas as gpd
import pyproj
from typing import Any, Iterable, Union, Tuple

import shapely.geometry

logger = logging.getLogger(__name__)


def bounds_to_left_top_right_bottom(bounds: Iterable[float], input_crs: Any):
    """
    bounds must be an iterable with four values (minx, miny, maxx, maxy)
    input_crs can be anything accepted by pyproj.CRS, such as: PROJ string, dictionary of PROJ params,
    authority string (epsg:4326), authority code as an int, a pyproj.CRS ...

    Returns the bounds as left, top, right, bottom, taking into account the axis_info in input_crs.
    """
    minx = bounds[0]
    miny = bounds[1]
    maxx = bounds[2]
    maxy = bounds[3]

    axis_info = pyproj.CRS(input_crs).axis_info

    if (axis_info[0].name == "Easting" and axis_info[0].direction == "east"
        and axis_info[1].name == "Northing" and axis_info[1].direction == "north") or \
            (axis_info[0].name == "Longitude" and axis_info[0].direction == "east"
             and axis_info[1].name == "Latitude" and axis_info[1].direction == "north"):
        # Most common situation (X/Y or lon/lat, Easting/Northing)
        left, top, right, bottom = minx, maxy, maxx, miny
    else:
        # TODO: Check which combination we have, and set left, top, right and bottom in accordance
        # Every possible combination between Easting/Northing, Northing/Easting,
        # Westing/Southing, Southing/Westing... I don't know if all of them are really possible
        raise NotImplementedError(f"This axis combination {axis_info} has not been implemented and tested yet")

    return left, top, right, bottom


def reproject_vector_file(input_file_path, output_file_path, dst_crs, layer=None, driver="GPKG"):
    """
    Reproject the vector file in input_file_path to output_file_path with a new CRS (dst_crs).
    By default it writes a GeoPackage.
    Simplistic, does not check anything...
    """
    data = gpd.read_file(input_file_path, layer=layer)
    data = data.to_crs(dst_crs)
    data.to_file(output_file_path, layer=layer, driver=driver)


def get_geodesic_area_perimeter_of_geom(geom: Union[shapely.geometry.Polygon, shapely.geometry.MultiPolygon],
                                        ellps:str ="WGS84") -> Tuple[float, float]:
    """
    Geodesic area and perimeter between the given points, in square meters and meters respectively.
    :param geom:
    :param ellps:
    :return:
    """
    geod = pyproj.Geod(ellps=ellps)
    area, perimeter = geod.geometry_area_perimeter(geom)
    return abs(area), abs(perimeter)


def get_geodesic_distance(lon1, lat1, lon2, lat2, ellps:str = "WGS84") -> float:
    """
    Geodesic distance between the given points, in meters.
    """
    geod = pyproj.Geod(ellps=ellps)
    fwazimuth, bwazimuth, dist = geod.inv(lon1, lat1, lon2, lat2)
    return dist