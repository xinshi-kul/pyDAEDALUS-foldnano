"""Turn a DAEDALUS routing (.cndo text) into helices and strands through foldnano.edit.

Layout of the result:

- A helix is one duplex of DAEDALUS's scaf_to_edge. Distances can't identify them: helices meet in
  line at some vertices, and base-pair spacing stretches on meshes with uneven edges. Its strand 0
  is the scaffold, which runs along +e3.
- An unpaired run (the 5-T vertex bulges) becomes a single-stranded extension past the end of the
  helix the strand is leaving, so every nucleotide and its order are kept.
- Crossovers in the routing leave the scaffold at about +107 degrees from e1 and the staples at
  about -107. Foldnano draws the strands 180 degrees apart, so strand 0 goes along e2 (+90), the
  closest fit to both.
"""

import math
from dataclasses import dataclass, field

from foldnano.region import BACKBONE_RISE_NM

ANGSTROM_PER_NM = 10.0
SCAFFOLD_COLOR = "#0066cc"
STAPLE_COLOR = "#f7931e"


@dataclass
class Nucleotide:
    id: int
    up: int
    down: int
    across: int
    seq: str


@dataclass
class Routing:
    nts: dict[int, Nucleotide]
    # Per base pair, in nm or as unit vectors.
    positions: dict[int, tuple[float, float, float]]
    triads: dict[int, tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]]
    pairs: dict[int, tuple[int, int]]
    nt_to_bp: dict[int, int] = field(default_factory=dict)

    def strands(self) -> list[list[int]]:
        """Nucleotide ids per strand, 5' to 3', longest first."""
        result = []
        for nt in self.nts.values():
            if nt.up != -1:
                continue
            strand = []
            current = nt.id
            while current != -1:
                strand.append(current)
                current = self.nts[current].down
            result.append(strand)
        if sum(map(len, result)) != len(self.nts):
            raise ValueError("Routing has a circular strand, which this plugin doesn't handle.")
        return sorted(result, key=len, reverse=True)


def parse_cndo(text: str) -> Routing:
    sections: dict[str, list[list[str]]] = {}
    current = None
    for line in text.splitlines():
        name = line.split(",")[0]
        if name in ("dnaTop", "dNode", "triad", "id_nt"):
            current = sections.setdefault(name, [])
        elif current is not None and line.strip():
            current.append(line.split(","))

    nts = {
        int(r[0]): Nucleotide(int(r[1]), int(r[2]), int(r[3]), int(r[4]), r[5])
        for r in sections["dnaTop"]
    }
    positions = {
        int(r[0]): tuple(float(v) / ANGSTROM_PER_NM for v in r[1:4]) for r in sections["dNode"]
    }
    triads = {
        int(r[0]): (
            tuple(map(float, r[1:4])),
            tuple(map(float, r[4:7])),
            tuple(map(float, r[7:10])),
        )
        for r in sections["triad"]
    }
    pairs = {int(r[0]): (int(r[1]), int(r[2])) for r in sections["id_nt"]}
    nt_to_bp = {nt: bp for bp, pair in pairs.items() for nt in pair}
    return Routing(nts, positions, triads, pairs, nt_to_bp)


def helices_from_route(routing: Routing, scaf_to_edge: list) -> list[list[int]]:
    """Base pairs per duplex, in scaffold 5' to 3' order, from DX_cage_design's scaf_to_edge."""
    # scaf_to_edge counts from 0, the .cndo from 1. A base pair has the id of its scaffold nucleotide.
    helices = [[int(nt) + 1 for nt in duplex] for edge in scaf_to_edge for duplex in edge]
    if sorted(bp for helix in helices for bp in helix) != sorted(routing.pairs):
        raise ValueError("scaf_to_edge and the .cndo don't list the same base pairs.")
    return helices


def _rotate_about(v: tuple[float, ...], axis: tuple[float, ...], degrees: float) -> tuple[float, float, float]:
    angle = math.radians(degrees)
    c, s = math.cos(angle), math.sin(angle)
    cross = (
        axis[1] * v[2] - axis[2] * v[1],
        axis[2] * v[0] - axis[0] * v[2],
        axis[0] * v[1] - axis[1] * v[0],
    )
    return tuple(c * v[i] + s * cross[i] for i in range(3))


def _bulge_room(routing: Routing, strands: list[list[int]]) -> int:
    longest = 0
    for strand in strands:
        run = 0
        for nt in strand:
            run = 0 if nt in routing.nt_to_bp else run + 1
            longest = max(longest, run)
    return longest


def routing_center(routing: Routing) -> tuple[float, float, float]:
    """Mean of the base-pair positions, in nm."""
    points = list(routing.positions.values())
    return tuple(sum(p[i] for p in points) / len(points) for i in range(3))


def build_design(
    part,
    region,
    routing: Routing,
    scaf_to_edge: list,
    center: tuple[float, float, float] = (0.0, 0.0, 0.0),
    twist_per_base_deg: float = 360 / 10.5,
) -> dict:
    """Add the routing to an empty region, centred on `center`. Returns counts for reporting."""
    strands = routing.strands()
    scaffold = strands[0]
    helices = helices_from_route(routing, scaf_to_edge)
    room = _bulge_room(routing, strands)
    longest = max(map(len, helices))

    needed = room + longest + room
    if region.max_base < needed:
        missing = needed - region.max_base
        region.increase_length_by(math.ceil(missing / region.step) * region.step)

    shift = tuple(center[i] - routing_center(routing)[i] for i in range(3))

    placement: dict[int, tuple[object, int]] = {}
    for helix in helices:
        mid_k = (len(helix) - 1) // 2
        mid_bp = helix[mid_k]
        e2, axis = routing.triads[mid_bp][1], routing.triads[mid_bp][2]
        mid_index = room + mid_k
        # Exact at the helix's middle: Foldnano's rise and twist differ a little from DAEDALUS's.
        origin = tuple(
            routing.positions[mid_bp][i] + shift[i] - axis[i] * mid_index * BACKBONE_RISE_NM
            for i in range(3)
        )
        vh = region.add_vh_at(origin, axis)
        radial0 = _rotate_about(e2, axis, -twist_per_base_deg * mid_index)
        vh.local_frame = {
            "origin": dict(zip("xyz", origin)),
            "axis": dict(zip("xyz", axis)),
            "radial0": dict(zip("xyz", radial0)),
        }
        for k, bp in enumerate(helix):
            placement[bp] = (vh, room + k)

    scaffold_ids = set(scaffold)
    for strand in strands:
        is_scaffold = strand is scaffold
        locations = _locate_strand(routing, strand, placement, 0 if is_scaffold else 1)
        segments = _segments_from(locations)
        oligo = part.create_oligo(
            segments,
            color=SCAFFOLD_COLOR if is_scaffold else STAPLE_COLOR,
            is_scaffold=is_scaffold,
        )
        seq_at = {(vh.index, bp): routing.nts[nt].seq for nt, (vh, bp) in zip(strand, locations)}
        for segment in oligo.segments:
            for bp in range(segment.low_index, segment.high_index + 1):
                segment.slot_at(bp).seq = seq_at[(segment.vh_index, bp)]

    return {
        "helices": len(helices),
        "strands": len(strands),
        "nucleotides": len(routing.nts),
        "base_pairs": len(routing.pairs),
        "max_base": region.max_base,
        "scaffold_length": len(scaffold_ids),
    }


def _locate_strand(routing: Routing, strand: list[int], placement, strand_index: int):
    """(vh, bp) for each nucleotide of a strand, with unpaired ones continuing past a helix end."""
    direction = 1 if strand_index == 0 else -1
    located: list = [None] * len(strand)
    last = None
    for i, nt in enumerate(strand):
        bp = routing.nt_to_bp.get(nt)
        if bp is not None:
            located[i] = last = placement[bp]
        elif last is not None:
            last = (last[0], last[1] + direction)
            located[i] = last
    # A strand that starts unpaired extends back from its first paired base.
    for i in range(len(strand) - 2, -1, -1):
        if located[i] is None:
            vh, bp = located[i + 1]
            located[i] = (vh, bp - direction)
    return located


def _segments_from(locations) -> list[tuple]:
    """Runs of neighbouring bases on one helix, as (vh, five, three)."""
    segments = []
    start = 0
    for i in range(1, len(locations) + 1):
        ends_run = (
            i == len(locations)
            or locations[i][0].index != locations[start][0].index
            or abs(locations[i][1] - locations[i - 1][1]) != 1
        )
        if ends_run:
            segments.append((locations[start][0], locations[start][1], locations[i - 1][1]))
            start = i
    return segments
