#!/usr/bin/env python3

import sys

sys.path.insert(1, "../../")
import test
import quantile as qt
from copy import deepcopy

TASK_ID = int(sys.argv[1])
gate_counter = 0
t_counter = 0


def mkgate(k, l, p):
    """
    Return a gate with name G_plaquete_p_unitary_i from qubit k (specified as a triple) to qubit l (specified as triple). Afer each time the function is called, i is increased by 1, starting at 0.
    The G_i are given in terms of standard gates in the paper. p labels the plaquette, either the first (0) or the second (2) plaquette
    """
    global gate_counter
    global t_counter
    name = f"G_plaquette_{p}_unitary_{gate_counter}"

    qk = qt.Qubit(*k)
    ql = qt.Qubit(*l)
    qubits = (qk, ql)

    gate = qt.Gate(name, qubits, t=t_counter)

    gate_counter += 1
    t_counter += 1
    return gate


def mk_basis_circ(decomposition):
    global gate_counter
    global t_counter
    if decomposition == "xyz":
        bes = [
            (0, 1),
            (2, 3),
            (1, 2),
            (0, 1),
            (1, 2),
            (2, 3),
            (1, 2),
            (0, 1),
            (1, 2),
            (0, 1),
            (2, 3),
            (1, 2),
            (0, 1),
            (2, 3),
            (0, 2),
            (1, 3),
            (0, 3),
            (1, 2),
        ]
    elif decomposition == "mcu":
        bes = [
            (0, 1),
            (2, 3),
            (0, 2),
            (1, 2),
            (2, 3),
            (1, 2),
            (2, 3),
            (1, 3),
            (0, 1),
            (0, 3),
            (0, 1),
            (0, 2),
            (0, 1),
            (0, 3),
            (0, 1),
            (0, 2),
            (0, 1),
            (2, 3),
        ]

    # Edges using correct seed notation
    es = [((0, 0, e[0]), (0, 0, e[1])) for e in bes]

    # Edges for the second plaquette in a basis circuit
    mapping = {
        (0, 1): ((0, 0, 2), (1, 0, 3)),
        (1, 2): ((0, -1, 3), (0, 0, 0)),
        (2, 3): ((0, 0, 0), (-1, 0, 1)),
        (3, 0): ((0, 1, 1), (0, 0, 2)),
        (0, 2): ((0, 0, 2), (1, 1, 0)),
        (3, 1): ((0, 0, 1), (1, -1, 3)),
    }
    _mapping = deepcopy(mapping)
    for key, val in _mapping.items():
        mapping[key[::-1]] = val[::-1]  # Just in case edges appear in the other order
    ses = [mapping[e] for e in bes]

    if decomposition == "xyz":
        # Convert the edges to gates with the correct names G_plaquette_p_unitary_i. See the paper for the definition of G_p_i in terms of elementary gates.
        gate_counter = 0
        gates = [mkgate(e[0], e[1], 0) for e in es[:-4]]
        gate_counter = 0
        t_counter = 0
        sgates = [mkgate(*e, 1) for e in ses[:-4]]
        _gates = []
        for i, gate in enumerate(gates):
            _gates += [gate, sgates[i]]
        gates = _gates

        # For the last 4 gates, the plaquettes go sequentially
        last_gates = [mkgate(e[0], e[1], 0) for e in es[-4:]]
        gate_counter = 14
        last_sgates = [mkgate(*e, 1) for e in ses[-4:]]
        gates += last_gates
        gates += last_sgates
    else:
        gate_counter = 0
        gates = [mkgate(e[0], e[1], 0) for e in es]
        gate_counter = 0
        sgates = [mkgate(*e, 1) for e in ses]
        gates += sgates

    bc = qt.BasisCirc(gates, name="Rokhsar-Kivelson_" + decomposition)

    for gate in bc.gates:
        print(gate)
    return bc


def route_rokhsar_kivelson(
    bc_size, bg_name, bg_size, merge_swaps, minimize_swaps, decomposition
):
    bc = mk_basis_circ(decomposition)
    bc = bc.resized(*bc_size)[0]

    bg = qt.BasisGraph.from_json(bg_name)
    bg = bg.resized(*bg_size)[0]

    t = qt.Transpiler(bc, bg)
    t.cyclic = True
    t.gate_dependencies = True
    t.merge_swaps = merge_swaps
    t.minimize_swaps = minimize_swaps

    sol = t.solve()
    return t


patches = [
    [(1, 1), "square", (2, 2)],
    [(2, 1), "square", (4, 2)],
    [(2, 2), "square", (4, 4)],
]

settings = []

for merge_swaps in [True, False]:
    for minimize_swaps in [False, True]:
        for decomposition in ["xyz", "mcu"]:
            for patch in patches:
                settings.append(
                    [
                        patch[0],
                        patch[1],
                        patch[2],
                        merge_swaps,
                        minimize_swaps,
                        decomposition,
                    ]
                )

args = settings[TASK_ID - 1]
print("Running test", TASK_ID - 1, "of", len(settings))
print("with settings", args)
t = route_rokhsar_kivelson(*args)
t.add_to_database("rokhsar-kivelson", fname="../solutions.pkl")
t.append_to_qasm_database("../solutions.qasm")
