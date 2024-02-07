import networkx as nx
import re

from typing import Sequence
import dggstools.auids.auids as auids


def is_a_subarea_rhealpix(auid1_comp_b64: str, auid2_comp_b64: str, nil_name: str = auids.DEFAULT_NIL_NAME) -> bool:
    """
    Takes two AUIDs and returns True if the first one is of an area which is spatially contained in the area
    defined by the second.
    If both AUIDs are generated from optimal sets of cuids, this always works (unless some bugs remain).
    If they are not optimal, it might give some false negatives.
    E.g.: lets assume N1, N2, N3 and N4 are the same area as N. If we generate an auid for N1, N2, N3 and N4, and
    another one for N, this function will return a false negative because the first list of cuids is not optimal
    (N1, N2, N3, N4 has an optimal version which is N).

    TODO: This is RHEALPix specific, for now at least; and it has just been lightly tested. It is a work in progress.
    TODO: Check optimality and maybe give a warning of potential false negatives if inputs are not optimal. Or make
    it explicit as a parameter to manifest if the inputs are optimal or not.
    """
    auid1_comp = auids.decode64(auid1_comp_b64)
    auid1 = auids.decompress_id_bytes(auid1_comp)
    auid2_comp = auids.decode64(auid2_comp_b64)
    auid2 = auids.decompress_id_bytes(auid2_comp)

    # Just checking if auid1 is a substring of auid2 would fail:
    # - The last substring of closing parenthesis must be removed from auid1 before checking. This is easy: just select
    #    everything up to the last nil_name: auid1[:auid1.rindex(nil_name)] in auid2
    # - CUIDs in rhealpix may start with N, S, P, O, Q or R. The AUID for [N1, S1] would fail to be seen as a subarea
    #   of the AUID for [N1, O1, S1], because the branch with the O1 in the prefix tree is "in the middle". A solution
    #  for rHEALPix (maybe for other DGGSs too if generalized) is considering separately the substrings starting
    #  with N, O, P, Q, R and S.

    first_char_substrings = {}  # A substring starting with each one of the allowed first chars
    substrings = re.split(r"([NOPQRS])", auid1)  # Parentheses make the split include the N,O,P,Q,R, or S in the result
    for idx in range(1, len(substrings), 2):  # First one is always the root char, we don't want it, so we start in 1
        first_char_substrings[substrings[idx]] = substrings[idx+1]

    result = True
    # For every substring (the one starting with N, the one with O, etc.) in the first auid, we check if it is
    # contained in the second one
    for fcs in first_char_substrings:
        fcs_substr = first_char_substrings[fcs]
        result = result and (fcs_substr[:fcs_substr.rfind(nil_name)] in auid2)

    return result


def _optimize_cuids_rhealpix_slow(cuids: Sequence[str], n_side: int) -> Sequence[str]:
    """
    This implementation is far slower than optimize_cuids_rehalpix. Kept just for testing (they must
    produce the same results for the same inputs; if not, at least one of them is wrong).
    """

    def substitute_full_levels(leafs):
        def all_are_leafs(nodes) -> bool:
            return_value = True
            for n in nodes:
                return_value = return_value and (n in T.predecessors(auids.NIL))
            return return_value

        nodes_to_remove = set()
        nodes_to_connect_to_nil = set()
        new_leaf_nodes = set()
        for leaf in leafs:
            # print(T.nodes[leaf]["source"])
            parent = next(T.predecessors(leaf))
            # T.successors(parent) is an iterator, if we want its len and then iterate over it again, copying
            # it to a list is the easiest way
            siblings = list(T.successors(parent))
            is_first_level = parent == auids.ROOT  # First level is N, O, P, Q, R or S. We don't want to remove them
            if not is_first_level and all_are_leafs(siblings) and len(siblings) == n_side * n_side:
                for s in siblings:
                    nodes_to_remove.add(s)
                nodes_to_connect_to_nil.add(parent)

        for n in nodes_to_connect_to_nil:
            T.add_edge(n, auids.NIL)
            new_leaf_nodes.add(n)

        for n in nodes_to_remove:
            T.remove_node(n)

        return new_leaf_nodes

    # Build the prefix tree
    T = auids.build_prefix_tree(sorted(cuids), True)
    # Just using a prefix tree removes the duplicates "for free"

    # And those cells that are occupying a full cell in the upper resolution now can be substituted by their parent.
    # (e.g. N1, N2, N3, and N4, as long as they are leafs and n_side==2 can be substituted by N).
    # For this:
    # If a "leaf" node (i.e., a node with an edge to the nil node), has 3 or 8 siblings (depending on the n_side),
    # that level is full and those nodes can be removed leaving only the parent node as new "leaf".
    # This must be repeated until no "leaf" nodes remain unchecked (including the ones we are adding).
    new_leafs = substitute_full_levels(T.predecessors(auids.NIL))
    while len(new_leafs) > 0:
        new_leafs = substitute_full_levels(new_leafs)

    # Return the optimized tree as a sequence of cuids
    return auids.prefix_tree_to_ids(T)


def optimize_cuids_rhealpix(cuids: Sequence[str], n_side: int) -> Sequence[str]:
    """
    Take a Sequence of cuids and returns an optimal version of it (optimal in the sense explained in the
    manuscript referred above). This is important because AUIDs are better (shorter, more "unique"...) when
    we know that they have been generated with optimal sets of cuids.

    This is RHEALPix specific, for now at least; and it has just been lightly tested. It is a work in progress.

    This becomes slow pretty fast as the resolution of the cuids increases. For a real administrative unit of some
    50.000 km² at resolution 8 (n_side 3), it takes about 30 s, while at resolution 9 is som 350 s.

    :param cuids:
    :return:
    """
    # Build the prefix tree
    T = auids.build_prefix_tree(sorted(cuids), True)
    # Just using a prefix tree removes the duplicates "for free"

    def _simplify(root_node: int):
        nodes_to_remove = set()
        nodes_to_connect_to_nil = set()
        nodes_to_check_again = set()
        for k, v in nx.bfs_successors(T, source=root_node):  # All descendents, not just neighbors, in breadth first order
            for n in v:
                # We don't need to check a node we are going to remove
                # And we don't need to check the "leafs" (predecessors of NIL) as they don't have
                # any successors to check
                if n not in nodes_to_remove and n not in T.predecessors(auids.NIL):
                    candidates_to_removal = set()
                    n_descendants = 0
                    highest_resolution_in_descendants = 0
                    current_resolution = nx.shortest_path_length(T, auids.ROOT, n) - 1
                    for kn, vn in nx.bfs_successors(T, source=n):
                        if len(vn) > 0:
                            if vn[0] != auids.NIL:
                                n_descendants += len(vn)
                                highest_resolution_in_descendants = max(highest_resolution_in_descendants,
                                                                        nx.shortest_path_length(T, auids.ROOT, vn[0]) - 1)
                                for candidate in vn:
                                    candidates_to_removal.add(candidate)

                    # Descendants of n is "complete" (max) if every level under it is complete
                    # (complete means that every possible cell up to the highest resolution in the descendants
                    # is present)
                    max_n_descendants = 0
                    for res in range(1, (highest_resolution_in_descendants - current_resolution) + 1):
                        max_n_descendants += (n_side * n_side) ** res

                    if n_descendants == max_n_descendants and max_n_descendants != 0:
                        for c in candidates_to_removal:
                            nodes_to_remove.add(c)
                        nodes_to_connect_to_nil.add(n)
                        # If I have removed a node, I need to check again the predecessors because
                        # I had discarded them as not complete, but after removal they might become
                        # complete themselves
                        for npred in T.predecessors(n):
                            nodes_to_check_again.add(npred)
        for n in nodes_to_remove:
            T.remove_node(n)
        for n in nodes_to_connect_to_nil:
            T.add_edge(n, auids.NIL)
        return nodes_to_check_again

    to_check_again = _simplify(auids.ROOT)
    while len(to_check_again) > 0:
        for n in to_check_again:
            # The original root is node 0, which is the predecessors of the resolution 0 cells (N, S etc.)
            # So the new root is not n, is its predecessor, unless n == ROOT
            if n == auids.ROOT:
                to_check_again = _simplify(n)
            else:
                to_check_again = _simplify(list(T.predecessors(n))[0])
    # Return the optimized tree as a sequence of cuids
    return auids.prefix_tree_to_ids(T)


def is_optimal_cuids_rhealpix(cuids: Sequence[str], n_side: int) -> bool:
    """
    Returns True if the cuids are optimal (in the sense explained in the manuscript referred above).
    This is generally pretty fast, even for large sequences of cuids. Much, much faster than optimizing them.
    (e.g., ~100 ms for this function for a sequence of cuids that requires ~130.000 ms to optimize in the
    same computer).
    The time for this function is very similar to the time that it takes to optimize an already optimal
    sequence, so it is not worth it to check optimality before attempting to optimize. Even when that optimization
    is unnecessary, you would not save much time (a few milliseconds).
    """

    T = auids.build_prefix_tree(sorted(cuids), True)

    # It will be optimal si thera are 0 "complete" nodes, being complete nodes
    # those who have all the possible children given the maximum resolution
    # reached by the tree under them
    for k, v in nx.bfs_successors(T, source=auids.ROOT):
        for n in v:
            if n not in T.predecessors(auids.NIL):
                n_descendants = 0
                highest_resolution_in_descendants = 0
                current_resolution = nx.shortest_path_length(T, auids.ROOT, n) - 1
                for kn, vn in nx.bfs_successors(T, source=n):
                    if len(vn) > 0:
                        if vn[0] != auids.NIL:
                            n_descendants += len(vn)
                            highest_resolution_in_descendants = max(highest_resolution_in_descendants,
                                                                    nx.shortest_path_length(T, auids.ROOT, vn[0]) - 1)
                max_n_descendants = 0
                for res in range(1, (highest_resolution_in_descendants - current_resolution) + 1):
                    max_n_descendants += (n_side * n_side) ** res

                if n_descendants == max_n_descendants and max_n_descendants != 0:
                    return False
    return True


class RHEALPixAUIDGenerator(auids.AUIDGenerator):
    """
    RHEALPix specific AUID generator. It will automatically optimize the cuids before generating the AUID.
    """
    def __init__(self, n_side: int,
                 pars: str = auids.DEFAULT_PARS, root_name: str = auids.DEFAULT_ROOT_NAME,
                 nil_name: str = auids.DEFAULT_NIL_NAME, with_opening_par: bool = auids.DEFAULT_WITH_OPENING_PAR,
                 with_trailing_pars: bool = auids.DEFAULT_WITH_TRAILING_PARS,
                 hash_digest_size: int = auids.DEFAULT_HASH_DIGEST_SIZE):

        super().__init__(pars, root_name, nil_name, with_opening_par, with_trailing_pars, hash_digest_size)
        self.n_side = n_side

    def generate_auid_hash_b64(self, cuids: Sequence[str]) -> (str, str):
        """
        Generates a b64 encoded string that can be safely used as an AUID for the given cuids.
        It also generates a b64 encoded hash string fot that AUID.

        If the given sequence is not optimal, it will be optimized before generating the AUID, so it generates
        the canonical AUID for the region defined by those cuids.
        """
        if not is_optimal_cuids_rhealpix(cuids, self.n_side):
            optimized_cuids = optimize_cuids_rhealpix(cuids, self.n_side)
        else:
            optimized_cuids = cuids
        return super().generate_auid_hash_b64(optimized_cuids)

    def cuids_from_auid_b64(self, auid_comp_b64: str) -> Sequence[str]:
        """
        Takes an AUID created with the generate_auid_hash_b64 method, and returns the list of cuids corresponding
        to that AUID. This list will be different from the original list of cuids used to generate the AUID if the
        original list was not optimal.
        """
        return super().cuids_from_auid_b64(auid_comp_b64)
