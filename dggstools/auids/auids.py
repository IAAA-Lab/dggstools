import logging
import uuid
import networkx as nx
import base64
import zlib
import hashlib


from networkx.classes import DiGraph
from typing import Tuple, Sequence

logger = logging.getLogger(__name__)

ROOT = 0  # In a nx prefix tree, the root node is always 0 and has "source" None.
NIL = -1  # In a nx prefix tree, the NIL node is the "pseudo-leaf" of the tree
# These are defaults that, I hope, do not collide with the symbols that may be used in different schemas for naming
# cells in DGGSs. You can use others if you want.
DEFAULT_ROOT_NAME = "^"
DEFAULT_NIL_NAME = "$"
DEFAULT_PARS = "¿?"
# These are the most efficient alternatives (for having shortest auids). In some cases we might want explicit opening
# parentheses or having the traling parentheses there (e.g., we want to produce a "closer-to-canonical" bp tree)
DEFAULT_WITH_OPENING_PAR = False
DEFAULT_WITH_TRAILING_PARS = False

DEFAULT_HASH_DIGEST_SIZE = 20


# Auxiliary functions to generate different string representations from string sets.

def build_prefix_tree(cuids: Sequence[str], are_sorted: bool = False) -> DiGraph:
    # A Prefix_Tree is essentially another fact_type for a trie
    if not are_sorted:
        sorted_cuids = sorted(cuids)
    else:
        sorted_cuids = cuids
    return nx.prefix_tree(sorted_cuids)


def generate_BP_from_tree(T: DiGraph, pars: str, root_name: str, nil_name: str,
                          with_opening_par: bool = DEFAULT_WITH_OPENING_PAR) -> str:
    """
    It returns a "balanced parenthesis" string representation of the tree T, including the source
    of the nodes (not only the parenthesis).
    :param pars: The style of parentheses can be chosen (a string with two chars, e.g. "()")
    :param root_name: The root node fact_type can be chosen (string with one char).
    :param nil_name: The fact_type for the nodes that are the string terminators in the tree/trie (string with one char)
    :param with_opening_par: You can choose if the BP has the opening parentheses or not (boolean). If you choose False
    the source parameter of the nodes of the Tree will be used instead (must be possible).
    """
    assert (len(pars) == 2)
    assert (len(root_name) == 1)
    assert (len(nil_name) == 1)
    # The NIL node makes a nx prefix tree technically not a tree, so this assertion would fail
    # assert is_tree(T)

    def _generate_BP_node(T: DiGraph, node: int, pars: str, nil_name: str, with_opening_par: bool) -> str:
        pref = pars[0] if with_opening_par else ""
        if node != NIL:
            candidate = pref + str(T.nodes[node]['source'])
            for son in T.adj[node]:
                candidate = candidate + _generate_BP_node(T, son, pars, nil_name, with_opening_par)
            candidate = candidate + pars[1]
        else:
            candidate = pref + nil_name + pars[1]
        return candidate

    pref = pars[0] if with_opening_par else ""
    s = pref + root_name
    for son in T.adj[ROOT]:
        s = s + _generate_BP_node(T, son, pars, nil_name, with_opening_par)
    s = s + pars[1]
    return s


# We define encode and decode functions just to makes sure that we are doing the same thing all the time

def encode64(bytes_obj: bytes) -> str:
    return base64.urlsafe_b64encode(bytes_obj).decode('utf-8')


def decode64(obj: str) -> bytes:
    return base64.urlsafe_b64decode(obj)


def hash_id(idstr: str, digest_size: int = DEFAULT_HASH_DIGEST_SIZE) -> (bytes, str):
    """
    Takes an id (a string) and returns a blake2b hash with digest_size, encoded as
    a utf-8 urlsafe base64 encoded string
    """
    m = hashlib.blake2b(idstr.encode(), digest_size=digest_size)
    m_bytes = m.digest()
    hash_b64_str = encode64(m_bytes)
    return m_bytes, hash_b64_str


def compress_id(idstr: str) -> bytes:
    """
    Returns a compression (zlib; as a bytes object) of id (must be a utf-8 string).
    """
    id_bytes = idstr.encode("utf-8")
    id_comp = zlib.compress(id_bytes)
    return id_comp


def decompress_id_bytes(id_bytes: bytes) -> str:
    """
    Returns a utf-8 string which is id_bytes (bytes object) decompressed.
    """
    id_decomp = zlib.decompress(id_bytes)
    return id_decomp.decode("utf-8")


# Auxiliary functions to generate a set of strings from its unique string representation.


def bp_auid_to_preffix_tree(bp_auid: str, pars: str, nil_name: str, with_opening_par: bool) -> DiGraph:
    # Init the tree with its root and the NIL "pseudo-leaf"
    t = nx.DiGraph()
    t.add_node(ROOT, source=None)
    t.add_node(NIL, source=NIL)
    # And populate with the contents of lbp
    current_node = ROOT
    populate_with_bp_auid(bp_auid, pars, nil_name, with_opening_par, t, ROOT, True)
    return t


def populate_with_bp_auid(bp_auid: str, pars: str, nil_name: str, with_opening_par: bool, t: DiGraph,
                          current_node: int, going_down: bool):
    while len(bp_auid) > 0:
        if with_opening_par and bp_auid[0] == pars[0]:
            going_down = True
            node_source, bp_auid = process_node_source(bp_auid[1:], pars, nil_name, '')
            new_node = str(uuid.uuid4())
            t.add_node(new_node, source=node_source)
            t.add_edge(current_node, new_node)
            current_node = new_node
        elif bp_auid[0] == pars[1]:
            if going_down:
                t.add_edge(current_node, NIL)
                going_down = False  # Add NIL only as a leaf
            current_node = next(t.predecessors(current_node))  # The only parent, this is a tree
            bp_auid = bp_auid[1:]
        else:
            going_down = True
            node_source = "" if bp_auid[0] == nil_name else bp_auid[0]
            bp_auid = bp_auid[1:]
            new_node = str(uuid.uuid4())
            t.add_node(new_node, source=node_source)
            t.add_edge(current_node, new_node)
            current_node = new_node


def process_node_source(bp_auid: str, pars: str, nil_name: str, node_source: str) -> (str, str):
    if bp_auid[0] != pars[0] and bp_auid[0] != bp_auid[1]:
        suf = "" if bp_auid[0] == nil_name else bp_auid[0]
        return process_node_source(bp_auid[1:], pars, nil_name, node_source + suf)
    else:
        return node_source, bp_auid


def prefix_tree_to_ids(T: DiGraph) -> Sequence[str]:
    ids = []
    for node in T.predecessors(NIL):
        idstr = ''
        while node != ROOT:
            idstr = str(T.nodes[node]['source']) + idstr
            node = next(T.predecessors(node))
        ids.append(idstr)
    return sorted(ids)


# The function generate_bp_auid_from_cuids is the one you need to take a number of cell unique identifiers (CUIDs),
# initialize a trie (preffix tree) with them, and get a balanced parentheses string representation of that trie,
# along with a hash and a compressed version of that string.

def generate_bp_auid_from_cuids(cuids: Sequence[str], pars: str = DEFAULT_PARS, root_name: str = DEFAULT_ROOT_NAME,
                                nil_name: str = DEFAULT_NIL_NAME, with_opening_par: bool = DEFAULT_WITH_OPENING_PAR,
                                with_trailing_pars: bool = DEFAULT_WITH_TRAILING_PARS) \
        -> Tuple[str, bytes, str, bytes, str, DiGraph]:
    """
    Takes a list of cell unique identifiers (cuids) and returns a tuple with:
    - A string representation of that list of cuids that can be used as an area unique identifier (AUID),
      based on a trie created from the sorted(cuids) and expressed as a balanced parenthesis string with or without the
      opening parenthesis
    - A blake2b hash of that string (as bytes, and as a b64 string)
    - A compressed version of that string (as bytes, and as a b64 string)
    - The trie (networkx object)

    If we call this function with two list of trivially equivalent (i.e.., they have repetitions) cuids
    such as ['N1', 'N1'] and ['N1'] this function will produce the same auid. This makes sense because both lists
    refer to the same area, but it is just a side effect of using prefix trees. Generally speaking, you
    should always optimize your cuids before passing them to this function. For a discussion about optimality of cuids,
    see <https://agile-online.org/images/conferences/2019/documents/short_papers/58_Upload_your_PDF_file.pdf>.

    TODO: Add a test of optimality of cuids, and an optimizer function.
    """
    # Sorting the ids is necessary to prevent that for instance ['N21', 'N22'] and ['N22', 'N21'] give
    # different results
    t = build_prefix_tree(sorted(cuids), True)
    auid = generate_BP_from_tree(t, pars, root_name, nil_name, with_opening_par)
    # Removing the trailing closing parentheses (i.e. those after the last NIL) is possible without
    # losing information (they can always be added if necessary, the information is contained in the remaining string).
    if not with_trailing_pars:
        auid = auid[:auid.rfind(nil_name)+1]

    hashed, hashed_b64 = hash_id(auid)
    auid_compressed = compress_id(auid)
    auid_compressed_b64 = encode64(auid_compressed)
    return auid, hashed, hashed_b64, auid_compressed, auid_compressed_b64, t


# The function generate_cuids_from_bp_auid is the one you need to take an area unique identifier (AUID) encoded
# as a balanced parenthesis string from a trie, and generate the cell unique identifiers (CUIDs) included in that area.


def generate_cuids_from_bp_auid(bp_auid: str, pars: str = DEFAULT_PARS, root_name: str = DEFAULT_ROOT_NAME,
                                nil_name: str = DEFAULT_NIL_NAME,
                                with_opening_par: bool = DEFAULT_WITH_OPENING_PAR,
                                with_trailing_pars: bool = DEFAULT_WITH_TRAILING_PARS) -> Sequence[str]:
    """
    Takes a string with a bp string auid and returns the list of cuids that are encoded
    in that auid.
    """
    if not with_trailing_pars:
        # We need to add the trailing parentheses before reconstructing the cuids
        # This is, essentially, counting one for each node in the tree and
        # subtracting one for each closing par. The remaining number is the trailing closing pars.
        count = 0
        for c in bp_auid:
            # This works both with and without opening parenthesis
            # (if there are opening pars, we need to count each node just once, so we don't count the opening par)
            if c not in pars:
                count += 1
            elif c == pars[1]:
                count -= 1
        bp_auid = bp_auid + pars[1] * count
    # The part with the root fact_type in the bp_auid is irrelevant: it will always be our
    # root node (with source None) in the tree, and it does not change at all the area encoded
    pref_len = 1 if with_opening_par else 0
    t = bp_auid_to_preffix_tree(bp_auid[len(root_name) + pref_len:-1], pars, nil_name, with_opening_par)
    return prefix_tree_to_ids(t)


class AUIDGenerator:
    """
    DGGS independent AUID generator.
    """
    def __init__(self, pars: str = DEFAULT_PARS, root_name: str = DEFAULT_ROOT_NAME, nil_name: str = DEFAULT_NIL_NAME,
                 with_opening_par: bool = DEFAULT_WITH_OPENING_PAR,
                 with_trailing_pars: bool = DEFAULT_WITH_TRAILING_PARS,
                 hash_digest_size: int = DEFAULT_HASH_DIGEST_SIZE):
        """
        Default values should be fine for any normal use case.
        """
        self.root_name = root_name
        self.nil_name = nil_name
        self.pars = pars
        self.with_opening_par = with_opening_par
        self.with_trailing_pars = with_trailing_pars
        self.hash_digest_size = hash_digest_size

    def generate_auid_hash_b64(self, cuids: Sequence[str]) -> (str, str):
        """
        Generates a b64 encoded string that can be safely used as an AUID for the given cuids.
        It also generates a b64 encoded hash string fot that AUID.
        The AUID will not be "canonical" for the region defined by the cuids, unless the given cuids are optimal for
        that region.
        """
        _, _, hashed_b64, _, auid_comp_b64, _ = generate_bp_auid_from_cuids(cuids, self.pars, self.root_name,
                                                                            self.nil_name, self.with_opening_par,
                                                                            self.with_trailing_pars)
        return auid_comp_b64, hashed_b64

    def cuids_from_auid_b64(self, auid_comp_b64: str) -> Sequence[str]:
        """
        Takes an AUID created with the generate_auid_hash_b64 method, and returns the list of cuids corresponding
        to that AUID. This list could be different from the original list of cuids used to generate the AUID if the
        original list was not optimal.
        :param auid_comp_b64:
        :return:
        """
        auid_comp = decode64(auid_comp_b64)
        auid = decompress_id_bytes(auid_comp)

        return generate_cuids_from_bp_auid(auid, self.pars, self.root_name, self.nil_name,
                                           self.with_opening_par, self.with_trailing_pars)

    def hash_b64_from_auid(self, auid_comp_b64: str) -> str:
        auid_comp = decode64(auid_comp_b64)
        auid = decompress_id_bytes(auid_comp)
        _, hash_b64_str = hash_id(auid, self.hash_digest_size)
        return hash_b64_str
