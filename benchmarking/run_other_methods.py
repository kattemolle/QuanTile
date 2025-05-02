#!/usr/bin/env python3

import networkx as nx
from copy import deepcopy
import sys

sys.path.insert(1, "../")
# sys.path.insert(1, "qiskit-ibm-transpiler")
import quantile as qt
import math
import cirq as cq
from time import time
import json
import qiskit as qk
from qiskit.transpiler import PassManager
import os

try:
    from qiskit_ibm_transpiler.ai.routing import AIRouting
except ModuleNotFoundError:
    print("AI routing not imported")

ROUTER = sys.argv[1]  # cirq or qiskit or qiskitAI
BASIS_CIRC = sys.argv[2]  # J1J2-line or J1J2-square


def load_J1J2_square():
    """
    Return J1J2 sqaure transpilerobject.
    """
    db = qt.load_solution_database("../circuits/solutions.pkl")
    tq = db["two-qudit"]
    flag = False
    for t in tq:
        if (
            t.basis_circ.name == "J1J2-square(2,2)"
            and t.basis_graph.name == "square(2,2)"
            and t.merge_swaps == False
            and t.cyclic == False
            and t.minimize_swaps == False
        ):
            assert flag == False
            # print("plotting", t)
            # t.plot_solution()
            transpiler = t
            flag = True

    return transpiler


def load_J1J2_line():
    """
    Return J1J2 transpiler object.
    """
    db = qt.load_solution_database("../circuits/solutions.pkl")
    tq = db["two-qudit"]
    flag = False
    for t in tq:
        if (
            t.basis_circ.name == "J1J2-line(4,1)"
            and t.basis_graph.name == "line(4,1)"
            and t.merge_swaps == False
            and t.cyclic == False
            and t.minimize_swaps == False
        ):
            assert flag == False
            # print("plotting", t)
            # t.plot_solution()
            transpiler = t
            flag = True

    return transpiler


def reconstruct_logical_circuit(c: qt.Circ, pqs: list[dict]):
    """
    Reconstuct logical circuit from a routed circuit c.
    """
    lqs = [{v: k for k, v in pqs[i].items()} for i in range(len(pqs))]
    gates = []
    for gate in c.gates:
        if gate.name != "rSWAP":
            _gate = deepcopy(gate)
            _gate.swap = False
            _qubits = []
            for _qubit in _gate.qubits:
                _qubits.append(lqs[gate.t][_qubit])
            _gate.qubits = tuple(_qubits)
            gates.append(_gate)

    lc = qt.Circ(gates)
    return lc


def to_router_circ(c):
    if ROUTER == "cirq":
        theta = 0.4 * math.pi
        phi = 2 * theta
        cirq_gates = []
        c.gates.sort(key=lambda x: x.t)
        for gate in c.gates:
            assert gate.name[0] == "G" and gate.swap == False
            qs = [cq.NamedQubit(str(qubit)) for qubit in gate.qubits]
            cirq_gate = cq.FSimGate(theta, phi).on(*qs)
            cirq_gates.append(cirq_gate)

        cirq_circuit = cq.Circuit(cirq_gates)
        return cirq_circuit
    elif ROUTER in ["qiskit", "qiskitAI"]:
        qubits = list(c.get_active_qubits())
        n = len(qubits)
        qubit_to_int = {qubits[i]: i for i in range(n)}
        qiskit_circ = qk.QuantumCircuit(n)
        c.gates.sort(key=lambda x: x.t)
        for gate in c.gates:
            assert gate.name[0] == "G" and gate.swap == False
            qs = gate.qubits
            qs = [qubit_to_int[q] for q in qs]
            qiskit_circ.cx(*qs)

        return qiskit_circ


def reconstruct_connectivity_graph(c):
    edges = []
    for gate in c.gates:
        qubits = [cq.NamedQubit(str(q)) for q in gate.qubits]
        edges.append(qubits)

    g = nx.Graph(edges)
    if ROUTER == "cirq":
        pass
    elif ROUTER in ["qiskit", "qiskitAI"]:
        g = list(nx.convert_node_labels_to_integers(g).edges())
        g += [edge[::-1] for edge in g]
        g = qk.transpiler.CouplingMap(g)

    return g


def count_swaps(c):
    if ROUTER == "cirq":
        count = 0
        for moment in c:
            for gate in moment:
                if isinstance(gate.gate, cq.SwapPowGate):
                    count += 1
                else:
                    assert isinstance(gate.gate, cq.FSimGate)
        return count
    elif ROUTER in ["qiskit", "qiskitAI"]:
        return c.count_ops()["swap"]


def get_active_qubits(c):
    qubits = set()
    for moment in c:
        for gate in c:
            for qubit in gate.qubits:
                qubits.add(qubit)

    return qubits


def router(transpiler, n, m):
    result = {}
    sol = transpiler.solution
    rbc = sol["routed_basis_circ"]
    pqs = sol["map_model"]
    routed_patch, pqs_routed_patch = rbc.get_completed_patch(
        n, m, phys_qubits=pqs, return_phys_qubits=True
    )
    connectivity_graph = reconstruct_connectivity_graph(routed_patch)
    qlc = reconstruct_logical_circuit(routed_patch, pqs_routed_patch)

    result["quantile_depth"] = sol["depth"]
    result["quantile_wall_clock"] = sol["wall_clock"]
    result["quantile_swaps"] = len(
        list(g for g in routed_patch.gates if g.name == "rSWAP")
    )
    result["quantile_active_qubits"] = len(routed_patch.get_active_qubits())

    clc = to_router_circ(qlc)
    start = time()
    if ROUTER == "cirq":
        router = cq.RouteCQC(connectivity_graph)
        routed_circuit = router(clc, lookahead_radius=32)
    elif ROUTER == "qiskit":
        routed_circuit = qk.transpile(
            clc,
            coupling_map=connectivity_graph,
            optimization_level=3,
            routing_method="sabre",
        )
    elif ROUTER == "qiskitAI":
        ai_passmanager = PassManager(
            [
                AIRouting(
                    coupling_map=connectivity_graph,
                    ai_optimization_level=3,
                    optimization_level=3,
                    local_mode=True,
                )
            ]
        )
        routed_circuit = ai_passmanager.run(clc)

    end = time()

    result["depth"] = len(routed_circuit)
    result["swaps"] = count_swaps(routed_circuit)
    result["active_qubits"] = len(get_active_qubits(routed_circuit))
    result["wall_clock"] = end - start

    return result


def mk_data(transpiler, linear, max_size, step_size):
    fname = f"{transpiler.name}_{ROUTER}.json"
    if os.path.isfile(fname):
        with open(fname, "r") as f:
            lines = json.load(f)
    else:
        lines = {}

    keys = [eval(key) for key in lines.keys()]
    max_key = max(keys)

    for i in range(max_key + 1, max_size + 1, step_size):
        print(i)
        if linear == True:
            j = 1
        else:
            j = i
        line = router(transpiler, i, j)
        print(line)
        lines[str(i)] = line

        with open(fname, "w") as f:
            json.dump(lines, f, indent=4)


if BASIS_CIRC == "J1J2-line":
    t = load_J1J2_line()
    mk_data(t, True, 250, 1)
elif BASIS_CIRC == "J1J2-square":
    t = load_J1J2_square()
    mk_data(t, False, 16, 1)
