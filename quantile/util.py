#!/usr/bin/env python3

from .base import Gate, Qubit, plot_dag, Circ, BasisGraph
from .transpiler import Transpiler
from copy import deepcopy
import networkx as nx


def reset_flags(c):
    """
    Reset all flags in DAG circuit c to False (in-place). Also needs to be run to create the flags (set to False) in the first place.
    """
    for n in c:
        if type(n) == Gate:
            for q in n.qubits:
                q.flag = False
        elif type(n) == tuple and type(n[0]) == Qubit:
            n[0].flag = False
        else:
            raise ValueError(n, type(n), type(n[0]))
    for e in c.edges(keys=True, data=True):
        e[3]["qubit"].flag = False


def update_path(c, node, old_q, new_q, final_map=None):
    """
    Given a DAG circuit c, return uc, where old_q in `node` is replaced by new_q, and likewise for all edges and nodes in the path descending from node. The flags of the visited qubits (in nodes and edges) are set to True.
    """
    on = deepcopy(node)  # Old node
    oes = [  # Outgoing edges of old node
        e
        for e in c.out_edges(on, keys=True, data=True)
        if e[3]["qubit"] == old_q and e[3]["qubit"].flag == False
    ]
    assert len(oes) in [0, 1]
    while len(oes) != 0:
        nn = deepcopy(on)  # New node
        new_qp = deepcopy(new_q)
        new_qp.flag = True
        if type(on) == Gate:
            nn.qubits = tuple(
                new_qp if q == old_q and q.flag == False else q for q in on.qubits
            )
        elif type(on) == tuple and type(on[0]) == Qubit and on[1] == "init":
            assert on[0] == old_q
            id = 0
            nn = (new_qp, "init", id)
            while nn in c.nodes:
                id += 1
                nn = (new_qp, "init", id)

        else:
            raise ValueError

        assert on in c, (on, c)
        cp = deepcopy(c)  # Seems to circumvent a bug in networkx
        c = nx.relabel_nodes(cp, {on: nn})

        assert (nn, oes[0][1], oes[0][2]) in c.edges
        nx.set_edge_attributes(
            c, {(nn, oes[0][1], oes[0][2]): {"qubit": deepcopy(new_qp)}}
        )

        # Prepare for next node
        on = oes[0][1]
        oes = [
            e
            for e in c.out_edges(on, keys=True, data=True)
            if e[3]["qubit"] == old_q and e[3]["qubit"].flag == False
        ]
        assert len(oes) in [0, 1]

    assert on[0] == old_q and on[1] == "final"
    new_qp = deepcopy(new_q)
    new_qp.flag = True
    id = 0
    nn = (new_qp, "final", id)
    while nn in c.nodes:
        id += 1
        nn = (new_qp, "final", id)
    assert on in c
    cp = deepcopy(c)  # Seems to circumvent a bug in networkx
    c = nx.relabel_nodes(cp, {on: nn})

    assert len(list(c.in_edges(nn))) == 1
    assert list(c.in_edges(nn, keys=True, data=True))[0][3]["qubit"] == new_qp

    return c


def check_init_final_map(c, im, fm):
    """
    Given DAG circuit c, return True if for every logical qubit lq, if you follow pq=im[lq] through the c, you end up at pq'=fm[lq].
    """

    def check_world_line(c, lq, im, fm):
        """
        Given a DAG circuit c, and a logical qubit lq, follow the logic qubit through the dag c. Return True if it ends up at the physical qubit consistent with fm.
        """
        lq = deepcopy(lq)
        lq.flag = False
        c = deepcopy(c)
        reset_flags(c)
        pq = im[lq]
        node = [
            node
            for node in c
            if type(node) == tuple and node[0] == pq and node[1] == "init"
        ]
        assert len(node) == 1
        node = node[0]
        oes = list(c.out_edges(node, data=True, keys=True))
        while len(oes) != 0:
            if type(node) == Gate and node.name == "rSWAP":
                node = [e[1] for e in oes if e[3]["qubit"] != pq]
                assert len(node) == 1
                node = node[0]
                pq = [e[3]["qubit"] for e in oes if e[3]["qubit"] != pq]
                assert len(pq) == 1
                pq = pq[0]
            else:
                node = [e[1] for e in oes if e[3]["qubit"] == pq]
                assert len(node) == 1
                node = node[0]
            oes = list(c.out_edges(node, data=True, keys=True))

        assert node[1] == "final"
        fpq = node[0]  # Final physical qubit
        return fpq == fm[lq]

    for lq in im:
        if not check_world_line(c, lq, im, fm):
            return False
    return True


def resolve_rswaps(rc):
    """
    Take the DAG routed circuit rc, remove all rSWAPS, and return the unrouted DAG circuit uc. Any rSWAP-node (inserted by routing) is removed and edges and nodes are relabeled accordingly. Merged swaps are not removed and assumed to be 'pulled out' to normal rSWAPs before rc is passed to this function.
    """
    uc = deepcopy(rc)
    clean = False
    while not clean:
        nodes = iter(list(uc.nodes))
        while True:
            try:
                node = next(nodes)
            except StopIteration:
                clean = True
                break
            if type(node) == Gate and node.name == "rSWAP":
                uc = resolve_rswap(node, uc)
                break

    return uc


def resolve_rswap(node, c):
    """
    Resolve the rSWAP `node` from DAG circuit c, return a deepcopy.
    """
    assert type(node) == Gate and node.name == "rSWAP"
    uc = deepcopy(c)
    ies = sorted(
        list(uc.in_edges(node, data=True, keys=True)),
        key=lambda e: e[3]["qubit"].to_tuple(),
    )  # In edges, sorted by qubit label
    oes = sorted(
        list(uc.out_edges(node, data=True, keys=True)),
        key=lambda e: e[3]["qubit"].to_tuple(),
    )  # Out edges, sorted by qubit label
    assert len(ies) == len(oes) == 2
    q0 = ies[0][3]["qubit"]
    q1 = ies[1][3]["qubit"]
    assert q0 == oes[0][3]["qubit"]
    assert q1 == oes[1][3]["qubit"]

    uc.remove_node(node)
    q0p = deepcopy(q0)
    q1p = deepcopy(q1)
    q0p.flag = True
    q1p.flag = True
    uc.add_edge(ies[0][0], oes[1][1], qubit=q0p)
    uc.add_edge(ies[1][0], oes[0][1], qubit=q1p)

    uc = update_path(uc, oes[1][1], q1, q0)
    other_target = [
        e for e in uc.out_edges(ies[1][0], keys=True, data=True) if e[3]["qubit"] == q1p
    ][0][1]
    assert other_target in uc
    uc = update_path(uc, other_target, q0, q1)

    reset_flags(uc)
    return uc


def dag_equal(c, cp, plot_dags=False, file1="c", file2="cp"):
    """
    Return true if the two dag circuits c, cp are equal. Only checks equality of all directed edges; ignores qubits on which no gates act and time attributes.
    """

    def clean_edges(c):
        """
        Return a cleaned up list of edges of DAG circuit c.
        """
        es = []
        for e in deepcopy(c.edges(data=True)):
            if type(e[0]) == tuple and type(e[1]) == tuple:
                assert e[0][1] == "init" and e[1][1] == "final"
                # Do not append the edge to es...
            else:
                ep = []
                for n in e[:2]:
                    if type(n) == Gate:
                        n.t = None
                        ep.append(n)
                    elif type(n) == tuple and type(n[0]) == Qubit:
                        ep.append(n[:2])
                    else:
                        raise ValueError
                ep.append({"qubit": e[2]["qubit"]})
                es.append(ep)
        return es

    _c = deepcopy(c)
    c = clean_edges(_c)
    _cp = deepcopy(cp)
    cp = clean_edges(_cp)

    if plot_dags == True:
        plot_dag(nx.MultiDiGraph(c), file1 + ".pdf")
        plot_dag(nx.MultiDiGraph(cp), file2 + ".pdf")

    cond = [c.count(el) == cp.count(el) for el in c]
    cond += [c.count(el) == cp.count(el) for el in cp]

    return all(cond)


def verification(t, n=4, m=4):
    """
    For a transpiler object t for which a solution was obtained, verify that the nxm patch generated from the routed circuit is equal to the nxm patch generated from the basis circuit by comparing their DAG representations.
    """
    assert t.solution != None, "First run Transiler.solve()"
    assert (
        t.solution["solved"] == True
    ), "There is no solution to be verified. Probably transpilation was impossible with the current transpiler settings."
    print("Verifying", t.name, "...", end=" ", flush=True)

    def unmerge_rswaps(rbc):
        rbc = deepcopy(rbc)
        for gate in rbc.gates:
            assert gate.t != None
            gate.t = gate.t * 2

        ess = []  # Extracted SWAPs
        for gate in rbc.gates:
            if gate.swap == True:
                swap = Gate("rSWAP", gate.qubits, t=gate.t + 1)
                ess.append(swap)
                gate.swap = False
        rbc.gates += ess
        return rbc

    bc = t.basis_circ
    c = bc.get_open_patch(n, m).to_dag()

    rbc = t.solution["routed_basis_circ"]
    phys_qubits = t.solution["map_model"]

    if t.merge_swaps == True:
        rbc = unmerge_rswaps(rbc)
        _phys_qubits = []
        for time_slice in phys_qubits:
            _phys_qubits.append(time_slice)
            _phys_qubits.append(time_slice)
        phys_qubits = _phys_qubits

    rc = rbc.get_completed_patch(
        n, m, phys_qubits=phys_qubits
    ).to_dag()  # This also test there are no collisions

    # Construct initial and final map for the patch
    map_model = t.solution[
        "map_model"
    ]  # (Map from logical qubits to physical qubit at every physical time T)
    im = map_model[0]  # Init map
    fm = map_model[-1]
    # For a n by m patch the init map and final must be enlarged
    eim = {}  # Enlarged init map
    efm = {}
    for lq, pq in im.items():
        for i in range(n):
            for j in range(m):
                eim[lq.translated(i, j)] = pq.translated(i, j)
    for lq, pq in fm.items():
        for i in range(n):
            for j in range(m):
                efm[lq.translated(i, j)] = pq.translated(i, j)

    if not check_init_final_map(rc, eim, efm):
        # qt.plot_dag(rc, "false_init_final_map.pdf")
        # print("Problem with init or final map")
        return False

    uc = resolve_rswaps(rc)

    # Update all qubit labels in the unrouted circuit according to the init map.
    reset_flags(uc)
    for lq, pq in eim.items():
        assert (pq, "init") in uc, (pq, "init")
        if lq != pq:
            uc = update_path(uc, (pq, "init"), pq, lq, fm)

    reset_flags(uc)

    if t.gate_dependencies == True:
        eq = dag_equal(c, uc, False, "pre-" + t.name, "post-" + t.name)
    else:
        c = Circ.from_dag(c)
        for gate in c.gates:
            gate.t = None
        uc = Circ.from_dag(uc)
        for gate in uc.gates:
            gate.t = None
            cond = [c.gates.count(gate) == uc.gates.count(gate) for gate in c.gates]
            cond += [c.gates.count(gate) == uc.gates.count(gate) for gate in uc.gates]
        eq = all(cond)

    return eq


def route_qsim(
    basis_circuit_name,
    basis_circuit_size,
    basis_graph_name,
    basis_graph_size,
    merge_swaps=True,
    cyclic=False,
    gate_dependencies=False,
    slice_depth=None,
    minimize_swaps=False,
):
    # Create basis circuit
    bg = BasisGraph.from_json(basis_circuit_name)
    bg = bg.resized(*basis_circuit_size)[0]
    bc = bg.to_basis_circ(long_gate_names=False)

    # Create basis graph
    bg = BasisGraph.from_json(basis_graph_name)
    bg = bg.resized(*basis_graph_size)[0]
    t = Transpiler(bc, bg)

    # Set transpiler options
    t.merge_swaps = merge_swaps
    t.cyclic = cyclic
    t.gate_dependencies = gate_dependencies
    t.slice_depth = slice_depth
    t.minimize_swaps = minimize_swaps

    t.solve()
    return t
