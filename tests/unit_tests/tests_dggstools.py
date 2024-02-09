import unittest
import logging
from dggstools.rhpx import rhpxutils
import dggstools.rhpx.utils.utils as utils
import rhealpixdggs.dggs as rhp


class DGGSToolsTestCase(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger("geo2dggs")
        # Configuration for the logs (including those of every library used that produces logs)
        logging.basicConfig(format="[%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s", level=logging.INFO)

    def test_closest_resolution(self):
        rdggs_helper = rhpxutils.RHEALPixDGGSHelper(
            rhp.RHEALPixDGGS(ellipsoid=rhp.WGS84_ELLIPSOID, north_square=1, south_square=0, N_side=3))
        idx, res = rdggs_helper.get_closest_resolution(128)
        self.assertAlmostEqual(res, 169.57389885298727)
        idx, res = rdggs_helper.get_closest_resolution(100)
        self.assertAlmostEqual(res, 56.52463295099575)

    def test_higher_resolution(self):
        rdggs_helper = rhpxutils.RHEALPixDGGSHelper(
            rhp.RHEALPixDGGS(ellipsoid=rhp.WGS84_ELLIPSOID, north_square=1, south_square=0, N_side=3))
        idx, res = rdggs_helper.get_closest_higher_resolution(128)
        self.assertAlmostEqual(res, 56.52463295099575)

    def test_lower_resolution(self):
        rdggs_helper = rhpxutils.RHEALPixDGGSHelper(
            rhp.RHEALPixDGGS(ellipsoid=rhp.WGS84_ELLIPSOID, north_square=1, south_square=0, N_side=3))
        idx, res = rdggs_helper.get_closest_lower_resolution(128)
        self.assertAlmostEqual(res, 169.57389885298727)

    def test_change_extension(self):
        name = "coco.jpg"
        self.assertEqual("coco.tif", utils.change_extension(name, "tif"))
        name2 = "coco.loco.jpgg"
        self.assertEqual("coco.loco.gtif", utils.change_extension(name2, "gtif"))

    def test_geodesic_size_of_raster_profile(self):
        class FakeCRS(object):
            def __init__(self):
                self.is_projected = False

        crs = FakeCRS()

        left = -9.64
        top = 44.24
        right= 3.53
        bottom = 35.66
        res_x = 0.0333333
        res_y = -0.033333
        width = (right - left) / res_x
        height = abs((top - bottom) / res_y)

        fake_profile = {
            "transform": [res_x, 0, left, 0, res_y, top],
            "width": width,
            "height": height,
            "crs": crs
        }

        geod_diag, res = rhpxutils.get_geodesic_size_from_raster_profile(fake_profile)
        # I measured the true values roughly on QGIS and then I logged the result here and the logged values
        # looked good enough, so I use them as true values for these assertions
        self.assertAlmostEqual(geod_diag, 1471245.1842697694)
        self.assertAlmostEqual(res, 3120.012427237816)

        left = -180
        top = 90
        right = 180
        bottom = -90
        res_x = 0.0333333
        res_y = -0.033333
        width = (right - left) / res_x
        height = abs((top - bottom) / res_y)

        fake_profile = {
            "transform": [res_x, 0, left, 0, res_y, top],
            "width": width,
            "height": height,
            "crs": crs
        }

        geod_diag, res = rhpxutils.get_geodesic_size_from_raster_profile(fake_profile)
        # I measured the true values roughly on QGIS and then I logged the result here and the logged values
        # looked good enough, so I use them as true values for these assertions
        self.assertAlmostEqual(geod_diag, 20003931.458625447)
        self.assertAlmostEqual(res, 1656.6676042015517)
    def test_parent_children_cellid(self):
        self.assertEqual("", rhpxutils.get_parent_cellid("N"))
        self.assertEqual("N", rhpxutils.get_parent_cellid("N1"))
        self.assertEqual("N2", rhpxutils.get_parent_cellid("N23"))
        self.assertEqual("O123", rhpxutils.get_parent_cellid("O1230"))

        rdggs_helper = rhpxutils.RHEALPixDGGSHelper(
            rhp.RHEALPixDGGS(ellipsoid=rhp.WGS84_ELLIPSOID, north_square=1, south_square=0, N_side=2))
        self.assertEqual(["N0", "N1", "N2", "N3"], rhpxutils.get_children_cellids("N", rdggs_helper.rdggs))
        self.assertEqual(["S10", "S11", "S12", "S13"], rhpxutils.get_children_cellids("S1", rdggs_helper.rdggs))

        rdggs_helper = rhpxutils.RHEALPixDGGSHelper(
            rhp.RHEALPixDGGS(ellipsoid=rhp.WGS84_ELLIPSOID, north_square=1, south_square=0, N_side=3))
        self.assertEqual(["N0", "N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8"], rhpxutils.get_children_cellids("N", rdggs_helper.rdggs))
        self.assertEqual(["S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18"],
                         rhpxutils.get_children_cellids("S1", rdggs_helper.rdggs))

    def test_ascendant_cellid_at_resolution(self):
        self.assertIsNone(rhpxutils.get_ascendant_cellid_at_resolution_idx("N", 0))
        self.assertEqual("N", rhpxutils.get_ascendant_cellid_at_resolution_idx("N1", 0))
        self.assertEqual("N", rhpxutils.get_ascendant_cellid_at_resolution_idx("N0821", 0))
        self.assertEqual("N0", rhpxutils.get_ascendant_cellid_at_resolution_idx("N0821", 1))

    def test_descendant_cellids_at_resolution(self):
        rdggs_helper = rhpxutils.RHEALPixDGGSHelper(
            rhp.RHEALPixDGGS(ellipsoid=rhp.WGS84_ELLIPSOID, north_square=1, south_square=0, N_side=2))
        self.assertEqual([], rhpxutils.get_descendant_cellids_at_resolution_idx("N", rdggs_helper.rdggs, 0))
        self.assertEqual(["N0", "N1", "N2", "N3"], rhpxutils.get_descendant_cellids_at_resolution_idx("N", rdggs_helper.rdggs, 1))
        self.assertEqual(["N00", "N01", "N02", "N03",
                               "N10", "N11", "N12", "N13",
                               "N20", "N21", "N22", "N23",
                               "N30", "N31", "N32", "N33",],
                         rhpxutils.get_descendant_cellids_at_resolution_idx("N", rdggs_helper.rdggs, 2))
        self.assertEqual(['N000', 'N001', 'N002', 'N003',
                          'N010', 'N011', 'N012', 'N013',
                          'N020', 'N021', 'N022', 'N023',
                          'N030', 'N031', 'N032', 'N033',
                          'N100', 'N101', 'N102', 'N103',
                          'N110', 'N111', 'N112', 'N113',
                          'N120', 'N121', 'N122', 'N123',
                          'N130', 'N131', 'N132', 'N133',
                          'N200', 'N201', 'N202', 'N203',
                          'N210', 'N211', 'N212', 'N213',
                          'N220', 'N221', 'N222', 'N223',
                          'N230', 'N231', 'N232', 'N233',
                          'N300', 'N301', 'N302', 'N303',
                          'N310', 'N311', 'N312', 'N313',
                          'N320', 'N321', 'N322', 'N323',
                          'N330', 'N331', 'N332', 'N333'],
                         rhpxutils.get_descendant_cellids_at_resolution_idx("N", rdggs_helper.rdggs, 3))
        self.assertEqual(["N00", "N01", "N02", "N03"],
                         rhpxutils.get_descendant_cellids_at_resolution_idx("N0", rdggs_helper.rdggs, 2))
        self.assertEqual(['N000', 'N001', 'N002', 'N003',
                          'N010', 'N011', 'N012', 'N013',
                          'N020', 'N021', 'N022', 'N023',
                          'N030', 'N031', 'N032', 'N033'],
                         rhpxutils.get_descendant_cellids_at_resolution_idx("N0", rdggs_helper.rdggs, 3))

    def test_descendant_cellids_up_to_resolution(self):
        rdggs_helper = rhpxutils.RHEALPixDGGSHelper(
            rhp.RHEALPixDGGS(ellipsoid=rhp.WGS84_ELLIPSOID, north_square=1, south_square=0, N_side=2))
        self.assertEqual([], rhpxutils.get_descendant_cellids_up_to_resolution_idx("N", rdggs_helper.rdggs, 0))
        self.assertEqual(["N0", "N1", "N2", "N3"], rhpxutils.get_descendant_cellids_up_to_resolution_idx("N", rdggs_helper.rdggs, 1))
        self.assertEqual(["N0", "N1", "N2", "N3",
                          "N00", "N01", "N02", "N03",
                          "N10", "N11", "N12", "N13",
                          "N20", "N21", "N22", "N23",
                          "N30", "N31", "N32", "N33", ],
                         rhpxutils.get_descendant_cellids_up_to_resolution_idx("N", rdggs_helper.rdggs, 2))

        self.assertEqual([ "N00", "N01", "N02", "N03"],
                         rhpxutils.get_descendant_cellids_up_to_resolution_idx("N0", rdggs_helper.rdggs, 2))
        self.assertEqual(rhpxutils.get_descendant_cellids_at_resolution_idx("N", rdggs_helper.rdggs, 1) +
                         rhpxutils.get_descendant_cellids_at_resolution_idx("N", rdggs_helper.rdggs, 2) +
                         rhpxutils.get_descendant_cellids_at_resolution_idx("N", rdggs_helper.rdggs, 3),
                            rhpxutils.get_descendant_cellids_up_to_resolution_idx("N", rdggs_helper.rdggs, 3))

    def test_ascendant_cellids_up_to_resolution(self):
        self.assertEqual([], rhpxutils.get_ascendant_cellids_up_to_resolution_idx("N", 0))
        self.assertEqual([], rhpxutils.get_ascendant_cellids_up_to_resolution_idx("N", 1))
        self.assertEqual(["N"], rhpxutils.get_ascendant_cellids_up_to_resolution_idx("N0", 0))
        self.assertEqual([], rhpxutils.get_ascendant_cellids_up_to_resolution_idx("N0", 1))
        self.assertEqual(["N", "N1", "N12", "N123", "N1231"], rhpxutils.get_ascendant_cellids_up_to_resolution_idx("N12313", 0))
        self.assertEqual(["N12", "N123", "N1231"],
                         rhpxutils.get_ascendant_cellids_up_to_resolution_idx("N12313", 2))

    def test_crs_to_rdggs(self):
        dggs = rhp.RHEALPixDGGS(ellipsoid=rhp.WGS84_ELLIPSOID, north_square=2, south_square=1, N_side=2)
        rdggs_helper = rhpxutils.RHEALPixDGGSHelper(dggs)
        rdggs_crs = rdggs_helper.rhealpixdef_to_pyproj_crs()
        dggs2 = rhpxutils.pyproj_crs_to_rdggs(rdggs_crs, 2)

        self.assertEqual(dggs.ellipsoid, dggs2.ellipsoid)
        self.assertEqual(dggs.south_square, dggs2.south_square)
        self.assertEqual(dggs.north_square, dggs2.north_square)
        self.assertEqual(dggs.N_side, dggs2.N_side)




if __name__ == '__main__':
    unittest.main()
