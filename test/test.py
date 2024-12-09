#!/usr/bin/env python3
"""
Run a variety of routings. Each time, test that a 5x5 patch of the routed circuit actually implements the 5x5 patch of the input circuit. This is done by transforming both the input circuit and routed circuit to a DAG circuit representation, resolving the rSWAPS in the routed circuit, and testing that the resulting DAG is equal to the DAG of a 5x patch of the input circuit.

Also tests all solutions in `circuits/solutions.pkl`. 

This verification step can be performed separately by calling `verification(t)` for some (solved) transpiler object `t`.

The file also contains some manual sub-tests that can be run separately for debugging.
"""
import sys

sys.path.insert(1, "../")
import quantile as qt
from quantile import util as ut
from copy import deepcopy
import numpy as np


def test_resolve_rswap():
    q0 = qt.Qubit(0, 0, 0)
    q1 = qt.Qubit(0, 0, 1)
    gate = qt.Gate("rSWAP", (q0, q1), 0)
    gate1 = qt.Gate("rSWAP", (q0, q1), 1)
    # gate2 = qt.Gate("rSWAP", (q1, q2), 2)
    c = qt.Circ([gate, gate1]).to_dag()
    qt.plot_dag(c, "one_rSWAP.pdf")

    c = ut.resolve_rswap(gate, c)
    qt.plot_dag(c, "0_rSWAP_resolved.pdf")

    gate1 = qt.Gate("rSWAP", (q1, q0), 1)

    c = ut.resolve_rswap(gate1, c)
    qt.plot_dag(c, "1_rSWAP_resolved.pdf")

    # gate2 = qt.Gate("rSWAP", (q0, q2), 2)

    # c = resolve_rswap(gate2, c)
    # qt.plot_dag(c, "2_rSWAP_resolved.pdf")


def test_resolve_rswaps():
    q0 = qt.Qubit(0, 0, 0)
    q1 = qt.Qubit(0, 0, 1)
    q2 = qt.Qubit(1, 0, 0)
    gate = qt.Gate("rSWAP", (q0, q1), 0)
    gate1 = qt.Gate("rSWAP", (q0, q1), 1)
    gate2 = qt.Gate("rSWAP", (q1, q2), 2)
    c = qt.Circ([gate, gate1, gate2]).to_dag()
    qt.plot_dag(c, "one_rSWAP.pdf")

    c = ut.resolve_rswaps(c)
    qt.plot_dag(c, "one_rSWAP_resolved.pdf")


def test_update_path():
    bg = qt.BasisGraph.from_json("line")
    bg = bg.resized(2, 1)[0]
    bc = bg.to_basis_circ("HEIS")
    c = bc.get_open_patch(4, 1).to_dag()
    qt.plot_dag(c, "line.png")

    old_q = qt.Qubit(2, 0, 1)
    new_q = qt.Qubit(9, 9, 9)
    ut.reset_flags(c)
    c = ut.update_path(c, (qt.Qubit(2, 0, 1), "init"), old_q, new_q)
    old_q = qt.Qubit(2, 0, 0)
    c = ut.update_path(c, (qt.Qubit(2, 0, 0), "init"), old_q, new_q)
    qt.plot_dag(c, "line_update_path.png")


def test_verification(t):
    """
    Test the tests. If even a single random gate or rSWAP is added to the routed basis circuit rbc, verification that rbc implements the basis circuit bc should return False.
    """
    print("testing test ...", end=" ", flush=True)

    def inject_random_gate(rbc, phys_qubits):
        rbc = deepcopy(rbc)

        ri = np.random.randint

        smax = rbc.get_smax()
        tmax = rbc.get_tmax()
        n = ri(1, 3)
        t = ri(tmax + 1)
        q0 = np.random.choice(
            list(phys_qubits[t].values())
        )  # Random physical qubit that holds a logical qubit.
        if n == 1:
            qubits = (q0,)
        else:
            q1 = q0
            while q1 == q0:
                q1 = qt.Qubit(ri(-1, 2), ri(-1, 2), ri(smax + 1))
            qubits = (q0, q1)
            if ri(2) == 1:
                qubits = qubits[::-1]

        if ri(2) == 1:
            name = "rSWAP"
        else:
            name = "gate"
        g = qt.Gate(name, qubits, t)

        for gate in rbc.gates:
            if gate.t >= t:
                gate.t += 1
        phys_qubits.insert(t, deepcopy(phys_qubits[t]))
        gates = rbc.gates
        gates.append(g)
        rbc = qt.BasisCirc(
            gates,
            unique_seeds_checking=False,
            collision_checking=False,
            wrapped_circ_checking=False,  # Is not 'turned on' in return statement because it does not return a Bool
        )

        return rbc, phys_qubits

    t = deepcopy(t)
    rbc = t.solution["routed_basis_circ"]
    phys_qubits = t.solution["map_model"]

    # The injection of the random gate should lead to an invalid circuit in the first place, raise an exception upon verification, or verification should fail.
    try:
        rbcp, phys_qubitsp = inject_random_gate(rbc, phys_qubits)
        t.solution["routed_basis_circ"] = rbcp
        t.solution["map_model"] = phys_qubitsp
        if not rbcp.unique_seeds() or not rbcp.no_collisions():
            return True
        else:
            v = qt.verification(t)
            if v is False:
                print(
                    "Routed circuit with random gate injected fails verification (as it should) without also throwing an exception"
                )
            return not v

    except:
        return True


def test_merged_swaps():
    q0 = qt.Qubit(0, 0, 0)
    q1 = qt.Qubit(0, 0, 1)
    q2 = qt.Qubit(1, 0, 0)
    gate = qt.Gate("HEIS", (q0, q1), 0)
    gate.swap = True
    gate1 = qt.Gate("HEIS", (q1, q2), 1)
    gate1.swap = True
    # gate2 = qt.Gate("rSWAP", (q1, q2), 2)
    c = qt.Circ([gate, gate1]).to_dag()
    qt.plot_dag(c, "merged_swap.pdf")

    uc = ut.resolve_rswaps(c)

    print(ut.dag_equal(c, uc))

    # c = resolve_rswap(gate, c)
    # qt.plot_dag(c, "0_rSWAP_resolved.pdf")

    # gate1 = qt.Gate("rSWAP", (q1, q0), 1)

    # c = resolve_rswap(gate1, c)
    # qt.plot_dag(c, "1_rSWAP_resolved.pdf")

    # gate2 = qt.Gate("rSWAP", (q0, q2), 2)

    # c = resolve_rswap(gate2, c)
    # qt.plot_dag(c, "2_rSWAP_resolved.pdf")


def test_export_qasm():
    q0 = qt.Qubit(0, 0, 0)
    q1 = qt.Qubit(0, 0, 1)
    q2 = qt.Qubit(1, 0, 0)
    gate0 = qt.Gate("h", (q0,), 0)
    gate1 = qt.Gate("cx", (q0, q1), 1)
    gate2 = qt.Gate("HEIS(0.13)", (q1, q2), 2)
    gate2.swap = True
    gates = [gate0, gate1, gate2]
    c = qt.Circ(gates)

    c.export_qasm()


def test_from_dag():
    import json

    with open("basis_graphs.json", "r") as f:
        db = json.load(f)

    for name in db:
        print(name)
        bg = qt.BasisGraph.from_json(name)
        bg = bg.resized(3, 3)[0]
        bc = bg.to_basis_circ("HEIS")
        dag = bc.to_dag()
        bc2 = qt.Circ.from_dag(dag)
        dag2 = bc2.to_dag()
        assert ut.dag_equal(dag, dag2)


def test():
    def route_append_plot_and_test(args):
        t = ut.route_qsim(*args)
        t.add_to_database("test")
        t.append_to_qasm_database("test_solutions_qasm.txt")
        # t.plot_solution("test")
        # print(t.solution["routed_basis_circ"].gates)
        assert qt.verification(t)
        assert test_verification(t)

    for merge_swaps in [True, False]:
        for cyclic in [False, True]:
            for gate_dependencies in [True, False]:
                for slice_depth in [None, 2, 4]:
                    for minimize_swaps in [False, True]:
                        route_append_plot_and_test(
                            [
                                "J1J2-line",
                                (3, 1),
                                "line",
                                (3, 1),
                                merge_swaps,
                                cyclic,
                                gate_dependencies,
                                slice_depth,
                                minimize_swaps,
                            ]
                        )

    for merge_swaps in [True, False]:
        for cyclic in [False, True]:
            for gate_dependencies in [True, False]:
                for slice_depth in [None, 2, 4]:
                    for minimize_swaps in [False, True]:
                        route_append_plot_and_test(
                            [
                                "kagome",
                                (1, 1),
                                "heavy-hex",
                                (1, 1),
                                merge_swaps,
                                cyclic,
                                gate_dependencies,
                                slice_depth,
                                minimize_swaps,
                            ]
                        )
                        route_append_plot_and_test(
                            [
                                "kagome",
                                (1, 1),
                                "square",
                                (2, 2),
                                merge_swaps,
                                cyclic,
                                gate_dependencies,
                                slice_depth,
                                minimize_swaps,
                            ]
                        )
                        route_append_plot_and_test(
                            [
                                "kagome",
                                (1, 1),
                                "square-octagon",
                                (1, 1),
                                merge_swaps,
                                cyclic,
                                gate_dependencies,
                                slice_depth,
                                minimize_swaps,
                            ]
                        )


def test_database():
    """
    Test all solutions from the solution database, and test that the test fails if the solution is corrupted.
    """
    db = qt.load_solution_database("circuits/solutions.pkl")
    for group in db:
        if group != "test":  # The group 'test' was already tested.
            for i, t in enumerate(db[group]):
                if t.solution["solved"] == True and t.solution["map_model"] != None:
                    v = qt.verification(t)
                    assert v
                    print(v)
                    v = test_verification(t)
                    assert v
                    print(v)
                    print()


if __name__ == "__main__":
    print("Routing and testing:")
    test()
    print()
    print("Testing database:")
    test_database()
