#!/usr/bin/env python3

import quantile as qt
from time import time


import sys

TASK_ID = int(sys.argv[1])


def mkedges(i, j, k):  # Make the edges of one rule 54 unitary acting on qubits i,j,k.
    return [(i, j), (k, i), (i, j), (k, i), (k, j)]


def mkgate(k, l, gate_counter, t_counter):
    """
    Return a gate with name G_i from qubit k (specified as a triple) to qubit l (specified as triple). Afer each time the function is called, i is increased by 1, starting at 0.
    The G_i are given in terms of standard gates in the paper.
    """
    name = "G_{}".format(gate_counter)

    qk = qt.Qubit(*k)
    ql = qt.Qubit(*l)
    qubits = (qk, ql)
    gate = qt.Gate(name, qubits, t=t_counter)

    return gate


def mk_basis_circ():
    he = [  # Hyperedges of the rule 54 unitary
        ((0, 0, 0), (0, 0, 1), (0, 0, 2)),
        ((0, 0, 2), (0, 0, 3), (1, 0, 0)),
        ((-1, 0, 3), (0, 0, 0), (0, 0, 1)),
        ((0, 0, 1), (0, 0, 2), (0, 0, 3)),
    ]

    edges = []
    for i, j, k in he:
        edges += mkedges(i, j, k)

    ### Make all gates
    gate_counter = 0
    t_counter = 0
    gates = []

    for e in edges:
        gates.append(mkgate(*e, gate_counter, t_counter))
        gate_counter += 1
        t_counter += 1

    return qt.BasisCirc(gates, "rule54")


def route_rule54(bc_size, bg_name, bg_size, merge_swaps, minimize_swaps):
    # Create basis circ
    bc = mk_basis_circ()
    bc = bc.resized(*bc_size)[0]

    # Create basis graph
    bg = qt.BasisGraph.from_json(bg_name)
    bg = bg.resized(*bg_size)[0]

    # Create transpiler and set options
    t = qt.Transpiler(bc, bg)
    t.cyclic = True
    t.gate_dependencies = True
    t.merge_swaps = merge_swaps
    t.minimize_swaps = minimize_swaps
    t.solve()
    return t


patches = [  # First size always for rule54.
    [(1, 1), "line", (4, 1)],
    [(2, 1), "line", (8, 1)],
    [(3, 1), "line", (12, 1)],
    [(4, 1), "line", (16, 1)],
    [(1, 1), "ladder", (2, 1)],
    [(2, 1), "ladder", (4, 1)],
    [(4, 1), "ladder", (8, 1)],
    [(1, 1), "square", (4, 1)],
    [(2, 1), "square", (4, 2)],
    [(2, 1), "square", (2, 4)],
    [(2, 1), "square", (8, 1)],
    [(3, 1), "square", (3, 4)],
    [(3, 1), "square", (4, 3)],
    [(3, 1), "square", (6, 2)],
    [(3, 1), "square", (2, 6)],
]

settings = []

for merge_swaps in [True, False]:
    for minimize_swaps in [False, True]:
        for patch in patches:
            settings.append(
                [
                    patch[0],
                    patch[1],
                    patch[2],
                    merge_swaps,
                    minimize_swaps,
                ]
            )

args = settings[TASK_ID - 1]
print("Running test", TASK_ID - 1, "of", len(settings))
print("with settings", args)
t = route_rule54(*args)
t.add_to_database("rule54", fname="../solutions.pkl")
t.append_to_qasm_database("../solutions.qasm")
