#!/usr/bin/env python3

"""
Run various routings for the the simulation of a general two-qudit hamiltonian.
"""
import sys

sys.path.insert(1, "../../")
import quantile as qt

TASK_ID = int(sys.argv[1])

patches = [
    ["J1J2-line", (16, 1), "ladder", (8, 1)],
    ["J1J2-square", (3, 3), "square", (3, 3)],
    ["J1J2J3-square", (3, 3), "square", (3, 3)],
    ["heavy-hex", (1, 1), "square", (3, 2)],
    ["heavy-hex", (1, 1), "square", (4, 2)],
    ["heavy-hex", (1, 1), "plus-square", (1, 1)],
    ["heavy-hex", (2, 1), "plus-square", (2, 1)],
    ["heavy-hex", (2, 1), "square", (5, 2)],
    ["shuriken", (1, 1), "square", (3, 2)],
    ["shuriken", (1, 1), "square", (3, 3)],
    ["shuriken", (2, 1), "square", (6, 2)],
    ["shuriken", (2, 1), "square", (4, 3)],
    ["triangular", (2, 2), "square", (2, 2)],
    ["triangular", (3, 3), "square", (3, 3)],
    ["triangular", (4, 3), "square", (4, 3)],
    ["honeycomb", (1, 1), "square", (2, 1)],
    ["square-octagon", (1, 1), "square", (2, 2)],
    ["square-octagon", (2, 1), "square", (4, 2)],
    ["square-octagon", (2, 2), "square", (4, 4)],
    ["kagome", (1, 1), "square", (2, 2)],
    ["kagome", (2, 1), "square", (3, 2)],
    ["kagome", (2, 1), "square", (2, 3)],
    ["kagome", (1, 2), "square", (3, 2)],
    ["kagome", (1, 2), "square", (2, 3)],
    ["kagome", (2, 2), "square", (6, 2)],
    ["kagome", (2, 2), "square", (4, 3)],
    ["cross", (1, 1), "square", (3, 4)],
    ["cross", (1, 1), "square", (4, 4)],
    ["star", (1, 1), "square", (3, 2)],
    ["star", (1, 1), "square", (3, 3)],
    ["star", (2, 1), "square", (4, 3)],
    ["ruby", (1, 1), "square", (3, 2)],
    ["ruby", (1, 1), "square", (3, 3)],
    ["ruby", (2, 1), "square", (4, 3)],
    ["trellis", (2, 1), "square", (2, 2)],
    ["trellis", (2, 2), "square", (2, 4)],
    ["trellis", (3, 2), "square", (3, 4)],
    ["snub-square", (1, 1), "square", (2, 2)],
    ["snub-square", (2, 1), "square", (4, 2)],
    ["snub-square", (2, 2), "square", (4, 4)],
    ["bridge", (1, 1), "square", (3, 2)],
    ["bridge", (2, 1), "square", (4, 3)],
    ["union-jack", (2, 2), "diamond-square", (2, 2)],
    ["union-jack", (2, 2), "square", (4, 2)],
    ["union-jack", (2, 2), "square", (3, 3)],
    ["asanoha", (2, 2), "square", (4, 3)],
    ["asanoha", (2, 2), "square", (4, 4)],
    ["kisrhombille", (1, 1), "square", (3, 2)],
    ["kisrhombille", (2, 1), "square", (4, 3)],
    ["tetrille", (1, 1), "square", (3, 2)],
    ["tetrille", (2, 1), "square", (4, 3)],
    ["dice", (1, 1), "square", (3, 1)],
    ["dice", (1, 1), "square", (2, 2)],
    ["dice", (1, 1), "diamond-square", (2, 1)],
    ["dice", (2, 1), "square", (3, 2)],
    ["floret-pentagonal", (1, 1), "square", (3, 3)],
    ["cairo-pentagonal", (1, 1), "square", (3, 2)],
    ["cairo-pentagonal", (2, 1), "square", (4, 3)],
    ["prismatic-pentagonal", (2, 1), "square", (3, 2)],
    ["prismatic-pentagonal", (2, 2), "square", (4, 3)],
]


settings = []
gate_dependencies = False
slice_depth = None

for merge_swaps in [True, False]:
    for cyclic in [False, True]:
        for minimize_swaps in [False, True]:
            for patch in patches:
                settings.append(
                    [
                        patch[0],
                        patch[1],
                        patch[2],
                        patch[3],
                        merge_swaps,
                        cyclic,
                        gate_dependencies,
                        slice_depth,
                        minimize_swaps,
                    ]
                )

args = settings[TASK_ID - 1]
print("Running test", TASK_ID - 1, "of", len(settings))
print("with settings", args)
t = qt.route_qsim(*args)
t.add_to_database("two-qudit", fname="../solutions.pkl")
t.append_to_qasm_database("../solutions.qasm")
# assert qt.verification(t)
