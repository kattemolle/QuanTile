#!/usr/bin/env python3
"""
Create and route the circuit for the quantum simulation of the Fermi-Hubbard model on a square grid.  Checked with manuscript on 27.08.24.
"""

import sys

sys.path.insert(1, "../../")
import quantile as qt

gate_counter = 0


def mkgate(k, l):
    """
    Return a gate with name U_i from qubit (0,0,k) to qubit (0,0,l). Afer each time the function is called, i is increased by 1, starting at 0.
    """
    global gate_counter
    name = "U_{}".format(gate_counter)

    qk = qt.Qubit(0, 0, k)
    ql = qt.Qubit(0, 0, l)
    qubits = (qk, ql)

    gate = qt.Gate(name, qubits, t=gate_counter)

    gate_counter += 1
    return gate


edges = [
    (0, 4),
    (4, 3),
    (4, 1),
    (0, 4),
    (4, 3),
    (4, 1),
    (0, 4),
    (2, 4),
    (4, 1),
    (4, 3),
    (2, 4),
    (4, 1),
    (4, 3),
    (2, 4),
]

gates = [mkgate(*edge) for edge in edges]

# Now add the last two gates, that are the only gates reaching outside the unit cell.

name = "U_{}".format(gate_counter)
t = gate_counter
q0 = qt.Qubit(0, 0, 2)
q1 = qt.Qubit(0, 1, 0)
qubits = (q0, q1)
gate = qt.Gate(name, (q0, q1), t=t)
gates += [gate]
gate_counter += 1

name = "U_{}".format(gate_counter)
q0 = qt.Qubit(0, 0, 1)
q1 = qt.Qubit(1, 0, 3)
qubits = (q0, q1)
gate = qt.Gate(
    name, qubits, t=t
)  # This gate and the previous can be applied simultaneaously.
gates += [gate]

bc = qt.BasisCirc(gates, name="fermi-hubbard")
bg = qt.BasisGraph.from_json("plus-square")

# These parameters do not need to be iterated over since no swaps need to be inserted...
t = qt.Transpiler(bc, bg)
t.cyclic = True
t.gate_dependencies = True
t.merge_swaps = False
t.minimize_swaps = True
sol = t.solve()
t.add_to_database("fermi-hubbard", fname="../solutions.pkl")
t.append_to_qasm_database("../solutions.qasm")

# Note that this gives the circuit for the even plaquettes once gates U_i are compiled to the hardware's native gates. To obtain the odd plaquettes, take the same circuit but change \alpha_T as described in the paper. Then, using Circ.translated, manually construct a 2x2 patch of the routed circuit, containing even and odd plaquettes and reseed it. This forms the basis circuit of the full solution.
