import unittest
import logging
from typing import Sequence

from dggstools.auids.auids import generate_bp_auid_from_cuids, generate_cuids_from_bp_auid, AUIDGenerator
from dggstools.auids.rhpx_auids import is_a_subarea_rhealpix, optimize_cuids_rhealpix, is_optimal_cuids_rhealpix,\
    RHEALPixAUIDGenerator


class AUIDsTestCase(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger("geo2dggs")
        # Configuration for the logs (including those of every library used that produces logs)
        logging.basicConfig(format="[%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s", level=logging.INFO)

    def test_encode_decode_cuids_without_repetition(self):
        def _test(cuids: Sequence[str], with_opening_par: bool):
            auid_bp, _, _, _, _, _ = generate_bp_auid_from_cuids(cuids, with_opening_par=with_opening_par)
            assert (generate_cuids_from_bp_auid(auid_bp, with_opening_par=with_opening_par) == sorted(cuids))

        _test(['N11', 'N12', 'N2', 'N3'], True)
        _test(['N11', 'N12', 'N2', 'N3'], False)
        _test(['N1', 'N21', 'N22', 'N23', 'N25', 'S12', 'S13'], True)
        _test(['N1', 'N21', 'N22', 'N23', 'N25', 'S12', 'S13'], False)
        _test(["[1,(0,0)]2", "[1,(0,1)]2", "[1,(1,1)]2"], True)
        _test(["[1,(0,0)]2", "[1,(0,1)]2", "[1,(1,1)]2"], False)
        _test(
            ["[1,(1,1)]1", "[1,(2,2)]2", "[1,(2,3)]2", "[1,(3,2)]2", "[1,(3,3)]2"], True)
        _test(
            ["[1,(1,1)]1", "[1,(2,2)]2", "[1,(2,3)]2", "[1,(3,2)]2", "[1,(3,3)]2"], False)
        _test(
            ['N(0,0)3', 'N(1,0)1', 'N(1,1)1', 'N(3,1)2', 'N(3,2)2', 'N(0,7)3', 'N(7,6)3',
             'N(12,5)4', 'N(3,5)4', 'N(6,5)4', 'N(1,15)4', 'N(7,7)4', 'N(6,12)4', 'N(9,8)4', 'N(0,9)4'], True)
        _test(
            ['N(0,0)3', 'N(1,0)1', 'N(1,1)1', 'N(3,1)2', 'N(3,2)2', 'N(0,7)3', 'N(7,6)3',
             'N(12,5)4', 'N(3,5)4', 'N(6,5)4', 'N(1,15)4', 'N(7,7)4', 'N(6,12)4', 'N(9,8)4', 'N(0,9)4'], False)

    def test_encode_decode_cuids_with_repetitions(self):
        def _test_pair(cuids: Sequence[str], equiv_cuids: Sequence[str], with_opening_par: bool):
            auid_bp, _, _, _, _, _ = generate_bp_auid_from_cuids(cuids, with_opening_par=with_opening_par)
            auid_bp2, _, _, _, _, _ = generate_bp_auid_from_cuids(equiv_cuids, with_opening_par=with_opening_par)
            assert (generate_cuids_from_bp_auid(auid_bp, with_opening_par=with_opening_par) ==
                    generate_cuids_from_bp_auid(auid_bp2, with_opening_par=with_opening_par))

        _test_pair(['N1', 'N1'], ['N1'], True)
        _test_pair(['N1', 'N1'], ['N1'], False)
        _test_pair(['N(1,0)1', 'N(1,0)1', 'N(1,1)1'], ['N(1,0)1', 'N(1,1)1'], True)
        _test_pair(['N(1,0)1', 'N(1,0)1', 'N(1,1)1'], ['N(1,0)1', 'N(1,1)1'], False)

    def test_cuids_order_not_relevant(self):
        def _test_pair(cuids1: Sequence[str], cuids2: Sequence[str]):
            bp_auid1, _, _, _, _, _ = generate_bp_auid_from_cuids(cuids1)
            bp_auid2, _, _, _, _, _ = generate_bp_auid_from_cuids(cuids2)
            assert (bp_auid1 == bp_auid2)

        _test_pair(['N22', 'N21'], ['N21', 'N22'])
        _test_pair(['[1,(1,0)]1', '[1,(0,1)]2', '[1,(0,1)]3'],
                   ['[1,(1,0)]1', '[1,(0,1)]3', '[1,(0,1)]2'])

    def test_AUID_generator(self):
        # This test illustrates the basic, and common, usage of the AUID generator that can be used for
        # any DGGS. That generator can't, for now, optimize sequences of cuids, so it can't guarantee
        # that it produces the canonical AUID for a given sequence of cuids.
        def _test(in_cuids: Sequence[str]):
            auid_comp_b64, hashed_b64 = a.generate_auid_hash_b64(in_cuids)
            cuids = a.cuids_from_auid_b64(auid_comp_b64)
            assert (sorted(in_cuids) == cuids)  # As long as there are no repetitions in in_cuids this must be true

        a = AUIDGenerator()
        _test(['N1'])
        _test(['N(1,0)1', 'N(1,0)2', 'N(1,1)1'])
        _test(["[1,(1,1)]1", "[1,(2,2)]2", "[1,(2,3)]2", "[1,(3,2)]2", "[1,(3,3)]2"])
        _test(['N(1,0)1', 'N(1,1)1', 'N(3,3)2', 'N(3,2)2', 'N(1,1)2', 'N(0,7)3', 'N(7,6)3', 'N(0,0)3',
               'N(12,5)4', 'N(3,5)4', 'N(6,5)4', 'N(1,15)4', 'N(7,7)4', 'N(6,12)4', 'N(9,8)4', 'N(0,9)4'])

    def test_RHEALPix_AUID_generator(self):
        # This test illustrates the basic, and common, usage of the AUID generator which is specific for
        # RHEALPix. This generator can optimize sequences of cuids, so it can guarantee that it produces
        # the canonical AUID for a given sequence of cuids.
        def _test(in_cuids: Sequence[str]):
            auid_comp_b64, hashed_b64 = a.generate_auid_hash_b64(in_cuids)
            cuids = a.cuids_from_auid_b64(auid_comp_b64)
            assert (optimize_cuids_rhealpix(in_cuids, a.n_side) == cuids)

        a = RHEALPixAUIDGenerator(n_side=2)

        _test(['N1'])
        _test(['N1', 'N21', 'N22', 'N23', 'N24', 'S12', 'S13'])
        _test(["N11", "N11", "N12", "N21", "N22", "N23", "N24", "N311", "N312", "N313", "N314",
               "N321", "N322", "N323", "N324", "N331", "N332", "N333", "N334", "N341", "N342", "N343"])

    def test_AUID_is_subarea_of_other(self):
        def _test(with_trailing_pars: bool):
            a = AUIDGenerator(with_trailing_pars=with_trailing_pars)
            auid1_comp_b64, _ = a.generate_auid_hash_b64(['N1', 'N21', 'N22', 'N23', 'N25', 'S12', 'S13'])
            auid2_comp_b64, _ = a.generate_auid_hash_b64(['N1', 'N21', 'N22', 'N23', 'N25', 'S12', 'S13', 'S14'])
            assert (is_a_subarea_rhealpix(auid1_comp_b64, auid2_comp_b64))

            auid1_comp_b64, _ = a.generate_auid_hash_b64(['N1', 'N21', 'N22', 'N23', 'N25', 'S12', 'S13', 'S14'])
            auid2_comp_b64, _ = a.generate_auid_hash_b64(
                ['N1', 'N21', 'N22', 'N23', 'N25', 'O143', 'O5135321', 'S12', 'S13', 'S14', 'S145', 'S2123', 'S2124'])

            assert (is_a_subarea_rhealpix(auid1_comp_b64, auid2_comp_b64))

            auid1_comp_b64, _ = a.generate_auid_hash_b64(
                ['N1', 'N21', 'N22', 'N23', 'N25', 'S12', 'S13', 'S14', 'O143'])
            auid2_comp_b64, _ = a.generate_auid_hash_b64(
                ['N1', 'N21', 'N22', 'N23', 'N25', 'O143', 'O5135321', 'S12', 'S13', 'S14', 'S145', 'S2123', 'S2124'])

            assert (is_a_subarea_rhealpix(auid1_comp_b64, auid2_comp_b64))

            auid1_comp_b64, _ = a.generate_auid_hash_b64(['N1', 'N21', 'N22', 'N23', 'N25', 'S12', 'S13', 'O143'])
            auid2_comp_b64, _ = a.generate_auid_hash_b64(
                ['N1', 'N21', 'N22', 'N23', 'N25', 'O143', 'O5135321', 'S12', 'S13', 'S14', 'S145', 'S2123', 'S2124'])

            assert (is_a_subarea_rhealpix(auid1_comp_b64, auid2_comp_b64))

            auid1_comp_b64, _ = a.generate_auid_hash_b64(['N1', 'N21', 'N22', 'N23', 'N25', 'S12', 'S13', 'S15'])
            auid2_comp_b64, _ = a.generate_auid_hash_b64(['N1', 'N21', 'N22', 'N23', 'N25', 'S12', 'S13', 'S14'])
            assert (not is_a_subarea_rhealpix(auid1_comp_b64, auid2_comp_b64))

            auid1_comp_b64, _ = a.generate_auid_hash_b64(['N1', 'N21', 'N22', 'N23', 'N25', 'S13', 'S14'])
            auid2_comp_b64, _ = a.generate_auid_hash_b64(
                ['N1', 'N21', 'N22', 'N23', 'N25', 'O143', 'O5135321', 'S12', 'S13', 'S14', 'S145', 'S2123', 'S2124'])
            assert (not is_a_subarea_rhealpix(auid1_comp_b64, auid2_comp_b64))

            auid1_comp_b64, _ = a.generate_auid_hash_b64(['N1', 'N21', 'N22', 'N23', 'N25', 'S12', 'S13', 'O1434'])
            auid2_comp_b64, _ = a.generate_auid_hash_b64(
                ['N1', 'N21', 'N22', 'N23', 'N25', 'O143', 'O5135321', 'S12', 'S13', 'S14', 'S145', 'S2123', 'S2124'])
            assert (not is_a_subarea_rhealpix(auid1_comp_b64, auid2_comp_b64))

        _test(True)
        _test(False)

    def test_AUID_optimizer(self):
        optimal = optimize_cuids_rhealpix(["N11", "N11", "N12", "N21", "N22", "N23", "N24"], n_side=2)
        assert (optimal == ["N11", "N12", "N2"])

        optimal = optimize_cuids_rhealpix(
            ["N11", "N11", "N12", "N21", "N22", "N23", "N24", "N311", "N312", "N313", "N314",
             "N321", "N322", "N323", "N324", "N331", "N332", "N333", "N334", "N341", "N342", "N343", "N344"], n_side=2)
        assert (optimal == ["N11", "N12", "N2", "N3"])

        optimal = optimize_cuids_rhealpix(
            ["N11", "N11", "N12", "N21", "N22", "N23", "N24", "N311", "N312", "N313", "N314",
             "N321", "N322", "N323", "N324", "N331", "N332", "N333", "N334", "N341", "N342", "N343"], n_side=2)
        assert (optimal == ["N11", "N12", "N2", "N31", "N32", "N33", "N341", "N342", "N343"])

        optimal = optimize_cuids_rhealpix(["N", "S", "O", "Q"], n_side=2)
        assert (optimal == ["N", "O", "Q", "S"])

        optimal = optimize_cuids_rhealpix(["N1", "N2", "N3", "N4", "S", "O", "Q"], n_side=2)
        assert (optimal == ["N", "O", "Q", "S"])

        optimal = optimize_cuids_rhealpix(["R11", "R12", "R13", "R141", "R142", "R143", "R144"], n_side=2)
        assert (optimal == ["R1"])

        optimal = optimize_cuids_rhealpix(["R11", "R12", "R13", "R141", "R142", "R143"], n_side=2)
        assert (optimal == ["R11", "R12", "R13", "R141", "R142", "R143"])

        optimal = optimize_cuids_rhealpix(["R11", "R12", "R13", "R14", "R15", "R16", "R17", "R18", "R19"], n_side=3)
        assert (optimal == ["R1"])

    def test_is_optimal_cuids(self):
        assert (not is_optimal_cuids_rhealpix(["N11", "N11", "N12", "N21", "N22", "N23", "N24"], n_side=2))
        assert (is_optimal_cuids_rhealpix(["N11", "N12", "N2"], n_side=2))

        assert (not is_optimal_cuids_rhealpix(["N11", "N11", "N12", "N21", "N22", "N23", "N24", "N311", "N312", "N313",
                                               "N314", "N321", "N322", "N323", "N324", "N331", "N332", "N333", "N334",
                                               "N341", "N342", "N343", "N344"], n_side=2))
        assert (is_optimal_cuids_rhealpix(["N11", "N12", "N2", "N3"], n_side=2))
        assert (not is_optimal_cuids_rhealpix(["N11", "N11", "N12", "N21", "N22", "N23", "N24", "N311", "N312", "N313",
                                               "N314", "N321", "N322", "N323", "N324", "N331", "N332", "N333", "N334",
                                               "N341", "N342", "N343"], n_side=2))
        assert (is_optimal_cuids_rhealpix(["N11", "N12", "N2", "N31", "N32", "N33", "N341", "N342", "N343"], n_side=2))
        assert (not is_optimal_cuids_rhealpix(["R11", "R12", "R13", "R14", "R15", "R16", "R17", "R18", "R19"],
                                              n_side=3))
        assert (is_optimal_cuids_rhealpix(("R1",), n_side=3))


if __name__ == '__main__':
    unittest.main()
