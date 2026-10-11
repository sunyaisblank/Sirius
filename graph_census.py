"""Read-only in-memory stage 7 / Endpoint graph census.

Reviewed source epoch: 7f028007a4935ffdefe5ee9c886b138299f50dd8.
Run from the Sirius checkout with python3 -B. No files or modules are emitted.
These are structural DAG counts and optimistic cuts, not machine-code identity,
admission proofs, safe instruction elisions, or measured performance.
Root owns preservation and disposal of this disposable diagnostic.
"""

from pathlib import Path
import hashlib
import struct

m = {"__name__": "retained_review_memory"}
source = Path("src/sirius/kernels/retained_program.py")
assert hashlib.sha256(source.read_bytes()).hexdigest() == (
    "3754f0911b2b1a3d58276b23f4468d61a6092c65a1b522d7b557e373c6e96516"
)
exec(compile(source.read_text(), str(source), "exec"), m)


def dependencies(index):
    op, a, b, c = m["ops"][index]
    return ([] if op in (0, 1) else [a]
            + ([] if op in (6, 7, 8) else [b])
            + ([c] if op in (10, 11) else []))


def closure(roots, stop=()):
    seen, todo, stop = set(), list(roots), set(stop)
    while todo:
        index = todo.pop()
        if index in seen:
            continue
        seen.add(index)
        if index not in stop:
            todo.extend(dependencies(index))
    return seen


def fingerprints():
    result = []
    for index, (op, a, b, c) in enumerate(m["ops"]):
        raw = (struct.pack("<II", op, a) if op in (0, 1) else
               struct.pack("<I", op)
               + b"".join(result[d] for d in dependencies(index)))
        result.append(hashlib.sha256(raw).digest())
    return result


def census(nodes):
    return {op: sum(m["ops"][i][0] == op for i in nodes)
            for op in range(12) if any(m["ops"][i][0] == op for i in nodes)}


for name, build, endpoint_build, profile in (
    ("general", "build_transport_program", "build_endpoint_program", "metric_profile"),
    ("Schwarzschild", "build_schwarzschild_transport_program",
     "build_schwarzschild_endpoint_program", "schwarzschild_metric_profile"),
):
    saved, transport_roots = m["compile_program"], []

    def capture(outputs, *args, **kwargs):
        transport_roots.extend(outputs)
        return saved(outputs, *args, **kwargs)

    m["compile_program"] = capture
    transport = m[build](True)
    m["compile_program"] = saved
    live = closure(transport_roots)
    row = [m["P"](m["cache"][(1, i, 0, 0)]) for i in range(45)]
    H, ell = m["chart_metric_profile"](row[4:8], row, row[44], m[profile])
    values = [H.v.i] + [value.v.i for value in ell]
    derivatives = [value.i for value in H.d] + [value.i for e in ell for value in e.d]
    profile_nodes = closure(values + derivatives)
    anchors = profile_nodes & live
    print(name, "profile", {
        "Transport_nodes": transport["instructions"],
        "value_roots_already_live": sum(i in live for i in values),
        "derivative_roots_already_live": sum(i in live for i in derivatives),
        "complete_profile_nodes": len(profile_nodes),
        "missing_profile_nodes": len(profile_nodes - live),
        "missing_profile_opcodes": census(profile_nodes - live),
        "existing_anchor_nodes": len(anchors),
        "existing_anchor_opcodes": census(anchors),
    })
    f = fingerprints()
    anchor_fingerprints = {f[i] for i in anchors}
    value_fingerprints = {f[i] for i in values if i in live}
    endpoint_roots = []

    def capture_endpoint(outputs, *args, **kwargs):
        endpoint_roots.extend(outputs)
        return saved(outputs, *args, **kwargs)

    m["compile_program"] = capture_endpoint
    endpoint = m[endpoint_build](True)
    m["compile_program"] = saved
    endpoint_live = closure(endpoint_roots)
    f = fingerprints()
    for label, matches in (("H_ell_values_only", value_fingerprints),
                           ("all_identical_geometry_anchors", anchor_fingerprints)):
        cut = {i for i in endpoint_live if f[i] in matches}
        removed = endpoint_live - closure(endpoint_roots, cut)
        replaced = removed | cut
        print(name, label, {
            "Endpoint_original_nodes": len(endpoint_live),
            "anchor_nodes": len(cut),
            "other_nodes_disappearing_if_cut": len(removed),
            "arithmetic_replaced": sum(m["ops"][i][0] not in (0, 1) for i in replaced),
            "divisions_replaced": sum(m["ops"][i][0] == 5 for i in replaced),
            "square_roots_replaced": sum(m["ops"][i][0] == 6 for i in replaced),
        })
