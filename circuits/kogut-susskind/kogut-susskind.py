#!/usr/bin/env python3

import sys

sys.path.insert(1, "../../")
import quantile as qt

import test


gate_counter = 0
t_counter = 0

# template
# ((,,),(,,)),

TASK_ID = int(sys.argv[1])


def kinetic_edges(a, b, c, d):
    """
    Output a list of edges representing gates of one edge of the kinetic hamiltonian along qubits a,b,c,d, with a placed on the qutrit.
    """
    kes = [(d, c)]  # kinetic edges

    okes = [  # other kinetic edges
        (c, b),
        (b, a),
        (c, b),
        (d, c),
    ]

    kes = kes + 4 * okes

    return kes


def all_kinetic_edges():
    """
    Output list of all edges in the basis circuit for the sim of the kinetic part of the Hamiltonian.
    """
    es = []
    terms = [
        ((0, 0, 1), (0, 0, 0), (0, 0, 6), (0, 0, 2)),  # Bottom left to bottom right
        ((0, 0, 4), (0, 1, 2), (0, 0, 6), (1, 0, 0)),  # Top left to top right
        ((0, 0, 5), (0, 0, 0), (0, 0, 6), (0, 1, 2)),  # Bottom left to top left
        ((0, 0, 3), (0, 0, 2), (0, 0, 6), (1, 0, 0)),  # Bottom right to top right
    ]
    for term in terms:
        es += kinetic_edges(*term)

    return es


def magnetic_edges(b, r, t, l):
    """
    Output a list of edges representing gates of one plaquette of the magnetic Hamiltonian along qubits b,r,t,l (bottom, right, top, left). This defines the edges prior to gate cancellation and merging.
    """
    mes = [
        (l, t),
        (t, r),
        (r, b),
        (t, b),
        (r, b),
        (l, b),
        (r, b),
        (t, b),
        (r, b),
        (l, b),
        (t, r),
        (l, t),
    ]  # magnetic edges
    mes = 8 * mes
    return mes


def all_magnetic_edges():
    """
    Output a list of all edges representing gates of both plaquettes of the magnetic Hamiltonian along qubits b,r,t,l (bottom, right, top, left). This is prior to gate cancellation and merging.
    """
    b, r, t, l = (0, 0, 1), (0, 0, 3), (0, 0, 4), (0, 0, 5)
    es = magnetic_edges(b, r, t, l)

    b, r, t, l = (-1, -1, 4), (0, -1, 5), (0, 0, 1), (-1, 0, 3)
    es += magnetic_edges(b, r, t, l)

    return es


def remove_gates(gates):
    """
    Remove redundant gates. Sometimes removing a gate means merging it into the previous two-qubit gate.
    """
    irlst = [11, 12, 35, 36, 59, 60, 83, 84, 13, 37, 61, 85, 24, 48, 72]
    _irlst = [i + 96 for i in irlst]
    irlst += _irlst
    rlist = ["G_{}^magnetic".format(i) for i in irlst]
    _gates = [gate for gate in gates if gate.name not in rlist]
    return _gates


def mkgate(part, k, l):
    """
    Return a gate with name G_i^part from qubit k (specified as a triple) to qubit l (specified as triple). Afer each time the function is called, i is increased by 1, starting at 0.
    The G_i^part are given in terms of standard gates in the paper.
    """
    global gate_counter
    global t_counter
    name = "G_{}^{}".format(gate_counter, part)

    qk = qt.Qubit(*k)
    ql = qt.Qubit(*l)
    qubits = (qk, ql)

    gate = qt.Gate(name, qubits, t=t_counter)

    gate_counter += 1
    t_counter += 1
    return gate


def mk_basis_circ():
    global gate_counter
    gates = []

    for e in all_kinetic_edges():
        gates.append(mkgate("kinetic", *e))

    gate_counter = 0
    for e in all_magnetic_edges():
        gates.append(mkgate("magnetic", *e))

    gates = remove_gates(gates)
    bc = qt.BasisCirc(gates, name="kogut-susskind")
    return bc


def route_kogut_susskind(
    bc_size, bg_name, bg_size, merge_swaps, minimize_swaps, slice_depth
):
    # Create basis circ
    bc = mk_basis_circ()
    for g in bc.gates:
        print(g)
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
    t.slice_depth = slice_depth

    t.solve()
    return t


# patches = [ # First line always kogut-susskind
#        [(2,1),'diamond-square',(4,2)],
#        [(2,1),'square',(4,4)],
#    ]

patches = [  # First line always kogut-susskind
    [(1, 1), "diamond-square", (2, 2)],
]

settings = []
for merge_swaps in [True, False]:
    for minimize_swaps in [True, False]:
        for slice_depth in [10, 20, 40]:
            for patch in patches:
                settings.append(
                    [
                        patch[0],
                        patch[1],
                        patch[2],
                        merge_swaps,
                        minimize_swaps,
                        slice_depth,
                    ]
                )

args = settings[TASK_ID - 1]
print("Running test", TASK_ID - 1, "of", len(settings))
print("with settings", args)
t = route_kogut_susskind(*args)
t.add_to_database("kogut-susskind", fname="../solutions.pkl")
t.append_to_qasm_database("../solutions.qasm")
