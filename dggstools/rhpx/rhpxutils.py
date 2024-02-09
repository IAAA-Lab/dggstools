from collections import namedtuple
from typing import List

import fiona.transform
import rasterio.crs
import shapely.geometry
import shapely.ops
from rasterio.transform import Affine

from rhealpixdggs.dggs import RHEALPixDGGS, Cell
from rhealpixdggs.ellipsoids import WGS84_A, WGS84_F, WGS84_ELLIPSOID
from dggstools.rhpx.utils.rasterutils import *

# ROBERT GIBB, ALEXANDER RAICHEV, AND MICHAEL SPETH. THE RHEALPIX DISCRETE GLOBAL GRID SYSTEM (2013)
RHEALPIX_MEAN_AREAL_DISTORTION = 1.178  # Min, max and median (which actually is 1.177) are equal to the mean, as
# rhealpix is an equiareal projection

RHEALPixDGGSNamedTuple = namedtuple("RHEALPixDGGSNamedTuple", ["ellipsoid", "n_side", "north_square", "south_square"])


def rdggs_to_namedtuple(rdggs: RHEALPixDGGS) -> RHEALPixDGGSNamedTuple:
    assert rdggs.ellipsoid.a == WGS84_A and rdggs.ellipsoid.f == WGS84_F, \
        "Only the WGS84 Ellipsoid can be used for now"
    return RHEALPixDGGSNamedTuple("WGS84", rdggs.N_side, rdggs.north_square,
                                  rdggs.south_square)

def namedtuple_to_rdggs(rdggs_namedtuple: RHEALPixDGGSNamedTuple) -> RHEALPixDGGS:
    assert rdggs_namedtuple.ellipsoid == "WGS84", "Only the WGS84 Ellipsoid can be used for now"
    return RHEALPixDGGS(ellipsoid=WGS84_ELLIPSOID, N_side=rdggs_namedtuple.n_side,
                        north_square=rdggs_namedtuple.north_square, south_square=rdggs_namedtuple.south_square)


def pyproj_crs_to_rdggs(crs: pyproj.CRS, N_side: int) -> RHEALPixDGGS:
    # n_side is not part of the rHEALPix projection, but it is part of the rHEALPIX DGGS. Here it must be a parameter
    crs_dict = crs.to_dict()
    assert crs_dict["proj"] == "rhealpix", "crs must be rHEALPix"
    assert crs.ellipsoid.name == "WGS 84" or crs.ellipsoid.name == "WGS84" or crs.ellipsoid.name == "WGS_84", \
        "Only the WGS84 Ellipsoid is supported"
    assert N_side == 2 or N_side == 3, f"N_side must be 2 or 3 but it is {N_side}"

    return RHEALPixDGGS(ellipsoid=WGS84_ELLIPSOID, N_side=N_side, north_square=crs_dict["north_square"],
                        south_square=crs_dict["south_square"])


def cellidstr_to_suid(cellid: str) -> List:
    return list(cellid[0]) + [int(digit) for digit in cellid[1:]]

def cellid_resolution_idx(cellid: str) -> int:
    return len(cellid) - 1

# TODO: TEST ALL THESE METHODS

def get_parent_cellid(cellid: str) -> str:
    if len(cellid) <= 1:
        return ""
    else:
        return cellid[:-1]

def get_ascendant_cellid_at_resolution_idx(cellid: str, res: int) -> str | None:
    cellid_res = cellid_resolution_idx(cellid)
    if cellid_res  <= res:
        return None
    else:
        if res < cellid_res:
            return cellid[:res+1]
        else:
            return None

def get_ascendant_cellids_up_to_resolution_idx(cellid: str, res: int) -> list[str]:
    cellid_res = cellid_resolution_idx(cellid)
    if cellid_res <= res:
        return []
    else:
        result = []
        for i in range(res, cellid_res):
            result.extend([cellid[:i+1]])
        return result


def get_children_cellids(cellid: str, rhpx: RHEALPixDGGS) -> List[str]:
    return [cellid + str(i) for i in range(rhpx.N_side**2)]

def get_descendant_cellids_at_resolution_idx(cellid: str, rhpx: RHEALPixDGGS, res: int) -> list[str]:
    """
    Returns the descendant cellids of cellid at exactly resolution res (not up to res, at res).
    """
    if res <= cellid_resolution_idx(cellid):
        return []
    elif res == cellid_resolution_idx(cellid) + 1:
        return get_children_cellids(cellid, rhpx)
    else:
        result = []
        for cid in get_descendant_cellids_at_resolution_idx(cellid, rhpx, res-1):
            result.extend(get_children_cellids(cid, rhpx))
        return result

def get_descendant_cellids_up_to_resolution_idx(cellid: str, rhpx: RHEALPixDGGS, res: int) -> list[str]:
    """
    Returns the descendant cellids of cellid up to resolution res (not at res, up to res). It does
    not include itself.
    """
    result = []
    for i in range(cellid_resolution_idx(cellid), res+1):
        result.extend(get_descendant_cellids_at_resolution_idx(cellid, rhpx, i))
    return result

def get_gdf_attrs_from_rhealpix_file(input_file_path: str) -> dict:
    result = {}
    profile = get_raster_profile(input_file_path)

    with rasterio.open(input_file_path) as raster:
        n_side = int(raster.tags()["n_side"])
        rdggs = pyproj_crs_to_rdggs(pyproj.CRS(profile["crs"]), n_side)
        rdggs_helper = RHEALPixDGGSHelper(rdggs)

        result["left"], result["top"], result["right"], result["bottom"], resx, resy = (
            get_bbox_from_raster_profile(profile))
        resolution_idx_x, _ = rdggs_helper.get_closest_resolution(abs(resx))
        resolution_idx_y, _ = rdggs_helper.get_closest_resolution(abs(resy))  # resy is often a negative number
        assert resolution_idx_x == resolution_idx_y, \
            f"{input_file_path} is not a proper rhealpix file. Its cells are not squares."

        result["resolution_idx"] = resolution_idx_x
        result["res"] = resx

        result["rhealpixdggs"] = {"n_side": rdggs.N_side,
                                  "north_square": rdggs.north_square,
                                  "south_square": rdggs.south_square,
                                  "max_areal_resolution": rdggs.max_areal_resolution,
                                  "max_resolution": rdggs.max_resolution,
                                  "ellipsoid": rdggs.ellipsoid}

        result["height"] = raster.height
        result["width"] = raster.width
        result["nbands"] = raster.count
        result["nodata"] = raster.nodata
        result["nodatavals"] = raster.nodatavals
        result["dtypes"] = raster.dtypes

    return result


def gdf_attrs_to_rdggs(gdf_attrs: dict) -> RHEALPixDGGS:
    return RHEALPixDGGS(ellipsoid=WGS84_ELLIPSOID,
                        N_side=gdf_attrs["rhealpixdggs"]["n_side"],
                        north_square=gdf_attrs["rhealpixdggs"]["north_square"],
                        south_square=gdf_attrs["rhealpixdggs"]["south_square"])


class RHEALPixDGGSHelper:

    def __init__(self, rdggs: RHEALPixDGGS):
        self.rdggs = rdggs

    def rhealpixdef_to_proj_string(self) -> str:
        assert self.rdggs.ellipsoid.a == WGS84_A and self.rdggs.ellipsoid.f == WGS84_F, \
            "Only the WGS84 Ellipsoid can be used for now"
        # The current value of ellipsoid.f in the library corresponds to that of the GRS80 ellipsoid
        # and not to the WGS84 one (at least according to the values in Proj.4). It does not seem to initialize
        # trouble right now (I have tried changing it and I have seen no changes) but maybe it could explain
        # some future problems... Proj.4 defaults to GRS80 for rhealpix.

        # The planar origin lon_0 and lat_0 is by default 0,0, so the planar coordinates 0,0 will fall on
        # the Q3 cell. Notice that this is different from some depictions of the planar rHEALPix in some
        # papers. We can change Lon_0 because proj supports that, but not lat_0
        assert self.rdggs.ellipsoid.lat_0 == 0, f"Proj does not support a lat_0 parameter " \
                                                f"(and {self.rdggs.ellipsoid.lat_0} was requested). This will not work."

        # by default is +ellps=GRS80, but that is OK as long as the WGS84 in the rhealpixdggs library
        # has the parameters of the GRS80. Besides that, I have not been able to detect differences
        # using one or the other
        if self.rdggs.ellipsoid.lon_0 == 0 and self.rdggs.ellipsoid.lat_0 == 0:
            return f"+proj=rhealpix +south_square={self.rdggs.south_square} +north_square={self.rdggs.north_square}"
        else:
            return f"+proj=rhealpix +south_square={self.rdggs.south_square} +north_square={self.rdggs.north_square} " \
                   f"+lon_0={self.rdggs.ellipsoid.lon_0}"

    def rhealpixdef_to_wkt_string(self) -> str:
        return pyproj.CRS(self.rhealpixdef_to_proj_string()).to_wkt()

    def rhealpixdef_to_proj_ellps(self) -> str:
        assert self.rdggs.ellipsoid.a == WGS84_A and self.rdggs.ellipsoid.f == WGS84_F, \
            "Only the WGS84 Ellipsoid can be used for now"
        return "WGS84"

    def rhealpixdef_to_pyproj_crs(self) -> pyproj.CRS:
        return pyproj.CRS(self.rhealpixdef_to_proj_string())


    def cell_widths_for_all_resolutions(self) -> List[float]:
        return [self.rdggs.cell_width(i) for i in range(self.rdggs.max_resolution)]

    def get_closest_higher_resolution(self, base_resolution: float) -> Tuple[int, float]:
        for i in range(self.rdggs.max_resolution):
            if self.rdggs.cell_width(i) < base_resolution:
                return i, self.rdggs.cell_width(i)

    def get_closest_lower_resolution(self, base_resolution: float) -> Tuple[int, float]:
        for i in range(self.rdggs.max_resolution):
            if self.rdggs.cell_width(i) < base_resolution:
                return i - 1, self.rdggs.cell_width(i - 1)

    def get_closest_resolution(self, base_resolution: float) -> Tuple[int, float]:
        for i in range(self.rdggs.max_resolution):
            if self.rdggs.cell_width(i) < base_resolution:
                higher = i, self.rdggs.cell_width(i)
                lower = i - 1, self.rdggs.cell_width(i - 1)
                if (lower[1] - base_resolution) < (base_resolution - higher[1]):
                    return lower
                else:
                    return higher

    def planar_boundary(self):
        # Generating a non planar boundary would not work (or at least is far more tricky).
        # The N and S cap corners are all aligned in geodetic coordinates, so that "squares" become lines...
        cell_geoms = []
        for cell in self.rdggs.grid(0):  # Every cell of resolution 0 (i.e., the faces of the cube)
            (xmin, xmax), (ymin, ymax) = cell.xy_range()
            cell_geoms.append(shapely.geometry.box(xmin, ymin, xmax, ymax))

        return shapely.ops.unary_union(cell_geoms)

    def project_and_clip_to_rhealpix(self, geom: dict, crs: str = None):
        """
        Project geom to rdggs, clip against the auids (the faces of the cube) and return the resulting geometry.

        geom is a fiona/shapely geometry
        crs is a proj4 string.
        """
        dst_crs = self.rhealpixdef_to_proj_string()
        unfolded_cube = self.planar_boundary()
        geom_rhealpix = fiona.transform.transform_geom(crs, dst_crs, geom)
        return shapely.geometry.shape(geom_rhealpix).intersection(unfolded_cube)

    def align_transform(self, transform: Affine, dst_resolution_idx: int) -> Affine:
        """
        Takes a transform not aligned to the rhealpix grid, and returns the closest one which is
        aligned.

        In theory, this could be substituted by something like this:
        transform, width, height = rasterio.warp.aligned_target(transform, width, height, dst_resolution)
        where width and height are taken from the output of calculate_default_transform.
        Although this solution seems to align the result properly in the horizontal axis, it does not seem to
        align properly in the vertical axis: the center of the output pixels do not correspond to the
        centroids of the cells (it seems to align them with the edges of those cells). And this in turn
        creates other problems interpreting each of those pixels as a cell (edges are problematic for this).
        So, until proven otherwise, we will assume that this apparently too complex solution is the right one.
        """
        # This code takes the top-left coordinates of transform, sees to which rhealpix cell they belong to, and uses
        # the closest vertex of that cell as the new top-left for the transform
        current_left = transform[2]
        current_top = transform[5]
        current_topleft_cell = self.rdggs.cell_from_point(dst_resolution_idx, (current_left, current_top))
        if current_topleft_cell is not None:
            new_left, new_top = self._get_closest_vertex_in_rhealpix_cell(current_left, current_top, current_topleft_cell)
        else:
            # if current_topleft_cell is None, it means that current_left, current_top falls outside the N square
            # (or maybe the S for some southern hemisphere datasets). We can realign taking as reference a cell that
            # falls on the Equator for the left:
            alternative_topleft_cell = self.rdggs.cell_from_point(dst_resolution_idx, (current_left, 0))
            new_left, _ = self._get_closest_vertex_in_rhealpix_cell(current_left, current_top, alternative_topleft_cell)

            # And with one that falls in the N square for the top
            c = Cell(self.rdggs, ['N', 0])
            (xmin, xmax), (ymin, ymax) = c.xy_range()
            # I need to add/subtract a small tolerance for this to work
            alternative_topleft_cell = self.rdggs.cell_from_point(dst_resolution_idx, (xmin+0.000001, current_top-0.000001))

            if alternative_topleft_cell is None:
                # Our final option is trying with the S square
                c = Cell(self.rdggs, ['S', 0])
                (xmin, xmax), (ymin, ymax) = c.xy_range()
                # I need to add/subtract a small tolerance for this to work

                alternative_topleft_cell = self.rdggs.cell_from_point(dst_resolution_idx,
                                                                      (xmin+0.000001, current_top-0.000001))
            _, new_top = self._get_closest_vertex_in_rhealpix_cell(current_left, current_top, alternative_topleft_cell)

        return Affine.translation(new_left - current_left, new_top - current_top) * transform

    def align_bounds(self, bounds:Tuple[float, float, float, float], dst_resolution_idx: int) -> Affine:
        """
        Takes a bounds tuple(left,top,right,bottom) not aligned to the rhealpix grid, and returns the
        closest one which is aligned using the align_transform method.
        """
        return self.align_transform(Affine(bounds[2]-bounds[0], 0, bounds[0], 0, bounds[3]-bounds[1], bounds[1]),
                                    dst_resolution_idx)


    def get_bbox_in_rhealpix(self, input_crs: rasterio.crs.CRS, left: float, top: float, right: float,
                             bottom: float) -> Tuple[float, float, float, float]:
        """
        Takes the left, top, right, bottom coordinates of a bounding box, coordinates in input_crs, and
        returns them reprojected to rhealpix as defined in rdggs (left, top, right, bottom).
        """
        dst_crs = self.rhealpixdef_to_proj_string()
        transformer = pyproj.Transformer.from_crs(input_crs, dst_crs, always_xy=True)
        left_rheal, top_rheal = transformer.transform(left, top)
        right_rheal, bottom_rheal = transformer.transform(right, bottom)

        # It is possible, I think, that the left,top and right,bottom corners in the original CRS are not
        # the left,top and right,bottom corners after transforming to rhealpix (because the N and S squares are
        # rotated). So we generate the other two corners and take the leftmost, rightmost, topmost and bottom-most values:
        left_right_candidate, top_bottom_candidate = transformer.transform(left, bottom)
        left_rheal = min(left_rheal, left_right_candidate)
        right_rheal = max(right_rheal, left_right_candidate)
        top_rheal = max(top_rheal, top_bottom_candidate)
        bottom_rheal = min(bottom_rheal, top_bottom_candidate)

        left_right_candidate, top_bottom_candidate = transformer.transform(right, top)
        left_rheal = min(left_rheal, left_right_candidate)
        right_rheal = max(right_rheal, left_right_candidate)
        top_rheal = max(top_rheal, top_bottom_candidate)
        bottom_rheal = min(bottom_rheal, top_bottom_candidate)

        return left_rheal, top_rheal, right_rheal, bottom_rheal

    def _get_closest_vertex_in_rhealpix_cell(self, left: float, top: float, cell: Cell) -> Tuple[float, float]:
        """
         Returns the vertex in cell which is closer (cartesian distance) to the point left, top.
         We are assuming plane cells.
         """

        def cartesian_dist(x1, y1, x2, y2):
            return ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5

        distances_and_vertices = [(cartesian_dist(left, top, x, y), x, y) for x, y in cell.vertices()]
        _, closest_x, closest_y = min(distances_and_vertices)
        return closest_x, closest_y



