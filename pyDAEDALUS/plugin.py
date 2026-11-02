"""Route a wireframe region's mesh with DAEDALUS and build the result through foldnano.edit."""

import glob
import math
import os
import pickle
import tempfile

import matplotlib

matplotlib.use("Agg")

from Automated_Design.DX_cage_design import DX_cage_design  # noqa: E402
from Automated_Design.ply_to_input import ply_to_input  # noqa: E402
from foldnano import PluginContext, ReferenceGeometry, ReferenceMesh  # noqa: E402
from foldnano.edit import load_part  # noqa: E402
from matplotlib import pyplot as plt  # noqa: E402

from convert_cndo import ANGSTROM_PER_NM, build_design, parse_cndo, routing_center  # noqa: E402

M13_LENGTH = 7249


def _route(mesh: ReferenceMesh, min_edge_bp: int) -> tuple[str, dict]:
    """The .cndo text and route info for a mesh. DAEDALUS only writes these to files."""
    with tempfile.TemporaryDirectory() as out:
        ply_path = os.path.join(out, "mesh.ply")
        with open(ply_path, "w") as f:
            f.write(ReferenceMesh(mesh.vertices, faces=mesh.faces).to_ply())
        coordinates, edges, faces, lengths, name, seq, xo = ply_to_input(
            ply_path, out, min_edge_bp, False
        )
        DX_cage_design(
            coordinates, edges, faces, lengths, name, seq, xo, [], [], False, out, 1,
            print_to_console=False,
        )
        plt.close("all")
        with open(glob.glob(os.path.join(out, "*.cndo"))[0]) as f:
            cndo = f.read()
        with open(glob.glob(os.path.join(out, "routeInfo_*.pickle"))[0], "rb") as f:
            route = pickle.load(f)
    return cndo, route


def _open_edges(mesh: ReferenceMesh) -> int:
    uses: dict[tuple[int, int], int] = {}
    for face in mesh.faces:
        for a, b in zip(face, (*face[1:], face[0])):
            key = (min(a, b), max(a, b))
            uses[key] = uses.get(key, 0) + 1
    return sum(count != 2 for count in uses.values())


def _diagonal(points) -> float:
    return math.dist(*(tuple(f(p[i] for p in points) for i in range(3)) for f in (min, max)))


def run(ctx: PluginContext) -> None:
    part = load_part(ctx.get_part_as_foldnano_json(), active_region_id=ctx.active_region_id)
    region = part.active_region
    if region is None or region.kind != "WIREFRAME":
        raise ValueError("Select a wireframe region first.")
    if region.virtual_helices:
        raise ValueError("This region already has helices. Use an empty region.")
    geometry = region.reference_geometry
    if geometry is None or geometry.mesh is None or not geometry.mesh.faces:
        raise ValueError("This region needs a reference mesh with faces.")
    mesh = geometry.mesh
    if _open_edges(mesh):
        raise ValueError("DAEDALUS needs a closed mesh: every edge must be shared by two faces.")

    routing_text, route = _route(mesh, int(ctx.get_param("min_edge_bp", 52)))
    routing = parse_cndo(routing_text)

    center = tuple(sum(v[i] for v in mesh.vertices) / len(mesh.vertices) for i in range(3))
    counts = build_design(part, region, routing, route["scaf_to_edge"], center)

    # DAEDALUS sets the size from the edge length, so the mesh follows the design to stay aligned.
    from_center = routing_center(routing)
    scaled = [
        tuple(float(v[i]) / ANGSTROM_PER_NM - from_center[i] + center[i] for i in range(3))
        for v in route["coordinates"]
    ]
    faces = [tuple(int(i) for i in face) for face in route["faces"]]
    region.set_reference_geometry(
        ReferenceGeometry(mesh=ReferenceMesh(scaled, faces=faces))
    )
    scale = _diagonal(scaled) / _diagonal([mesh.vertices[i] for face in mesh.faces for i in face])

    if counts["scaffold_length"] > M13_LENGTH:
        ctx.log(
            f"The scaffold is {counts['scaffold_length']} nt, longer than M13 ({M13_LENGTH} nt), "
            "so DAEDALUS made up a random sequence. It isn't a real scaffold."
        )
    ctx.log(
        f"{counts['helices']} helices, {counts['strands']} strands, "
        f"{counts['nucleotides']} nucleotides. The mesh was rescaled by {scale:.3g} to fit."
    )
    ctx.set_result(counts)
    ctx.set_result_mutation_log(part.mutation_log)
