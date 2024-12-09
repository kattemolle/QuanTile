#!/usr/bin/env python
"""
Library of constraints that can be added to the z3 solver.
"""
import itertools as it
import z3
from .base import Qubit


def mapping(
    solver, basis_graph, logic_qubits, phys_qubits, depth, Smax, init_map, final_map
):
    if init_map == None:
        # At T=0, logical qubits are fixed to their cell.
        for lq, pq in phys_qubits[0].items():
            solver.add(pq.x == lq.x, pq.y == lq.y)
            solver.add(0 <= pq.s, pq.s <= Smax)

    else:  # Use specified init map if one is specified
        assert type(init_map) == dict

        # assert set(init_map.keys()) == set(phys_qubits[0].keys())
        for lq, pq in phys_qubits[0].items():
            pgval = init_map[lq]
            solver.add(
                pq.x == pgval.x, pq.y == pgval.y, pq.s == pgval.s
            )  # This assumes that the init map is valid

    # At any later time, a logical qubit can only be in the base cell or in neighboring cells.
    r = len(phys_qubits)
    for T in range(1, r):  #     for T in range(r):
        for lq, pq in phys_qubits[T].items():
            solver.add(-1 <= pq.x, pq.x <= 1)
            solver.add(-1 <= pq.y, pq.y <= 1)
            solver.add(0 <= pq.s, pq.s <= Smax)

    # At any time, the map from logical qubits to physical qubits is injective.
    for T in range(r):
        for phys_qubit, phys_qubitp in it.combinations(phys_qubits[T].values(), 2):
            solver.add(z3.Not(phys_qubit == phys_qubitp))

    if final_map != None:
        assert type(final_map) == dict
        assert set(final_map.keys()) == set(phys_qubits[-1].keys())
        for lq, pq in phys_qubits[-1].items():
            pgval = final_map[lq]
            solver.add(pq.x == pgval.x, pq.y == pgval.y, pq.s == pgval.s)


def swap_effect(solver, basis_graph, logic_qubits, phys_qubits, rswaps, init_map, Smax):
    def nbrs(v, T):
        return {rswap for rswap in rswaps[T] if v in rswap.qubits}

    # If no swaps act on a physical qubit, at the next time step, the physical qubit holds the same logical qubit.
    for T, rswapsT in enumerate(rswaps):
        for lq in logic_qubits:
            for S in range(Smax + 1):
                for i in range(-1, 2):
                    for j in range(-1, 2):
                        v = Qubit(i, j, S)
                        qav = (
                            phys_qubits[T][lq] == v
                        )  # If the logical qubit is at vertex v in the mobility zone
                        andlist = [rswap.on == False for rswap in nbrs(v, T)]
                        no_swaps = z3.And(*andlist)  # and no swaps act on that vertex
                        qsav = (
                            phys_qubits[T + 1][lq] == v
                        )  # then After T, the logical qubit is still at v. Note phys_qubits must always have one more timestep than rswaps.
                        solver.add(z3.Implies(z3.And(qav, no_swaps), qsav))

    # If a swap acts on two physical qubits, the logical qubits they hold are interchanged.
    for T, rswapsT in enumerate(rswaps):
        for s, rswap in enumerate(rswapsT):
            for i, v in enumerate(rswap.qubits):
                for lq in logic_qubits:
                    qav = phys_qubits[T][lq] == v  # If logic qubit is at vertex v
                    swap_on = rswap.on  # and a swap acts on v
                    vp = rswap.qubits[(i + 1) % 2]
                    qavp = (
                        phys_qubits[T + 1][lq] == vp
                    )  # then at T+1, the logic qubit is at the over vertex vp of the edge the swap acts on.
                    solver.add(z3.Implies(z3.And(qav, swap_on), qavp))


def no_gate_gate_collisions(solver, phys_gates):
    """
    Two gates cannot act on the same qubit at the same time. (This is superfluous only if gate_dependencies are added to the solver).
    """
    # No seed collisions for any qubit in any pair of gates (including collisions of a gate with itself after wrapping).
    for lg, lgp in it.combinations(phys_gates.keys(), 2):
        for q0 in phys_gates[lg].qubits:
            for q1 in phys_gates[lgp].qubits:
                eqT = (
                    phys_gates[lg].t == phys_gates[lgp].t
                )  # If gates act at equal time
                eqs = z3.Not(
                    q0.s == q1.s
                )  # Then the seeds of q0 (from one gate) and q1 (from the other or the same gate) cannot be the same.
                solver.add(z3.Implies(eqT, eqs))


def no_gate_swap_collisions(solver, phys_gates, all_rswaps, merge_swaps, rswaps):
    """
    Iff `merge_swaps==True`, assume that two-qubit gates can be merged with swap gates. This will generally cause swap gates to be inserted at the exact space time locations of some two qubit gates. The actual merging has to be done in post-processing.
    """
    if merge_swaps == False:
        for T, rswapsT in enumerate(rswaps):
            for pg in phys_gates.values():
                to = pg.t == T  # Time overlap
                for swap in rswapsT:
                    for qubit in pg.qubits:
                        for v in swap.qubits:
                            so = qubit.s == v.s  # Space overlap after wrapping
                            nc = z3.Not(z3.And(to, so, swap.on))  # No collision
                            solver.add(nc)
    else:
        for T, all_rswapsT in enumerate(all_rswaps):
            for pg in phys_gates.values():
                to = pg.t == T  # Time overlap
                or_lst = [False]
                and_lst = []
                for swap in all_rswapsT:
                    if (
                        pg.n == 2
                    ):  # TODO this assumes a SWAP gate cannot be merged with two single-qubit gates.
                        eo = z3.Or(
                            z3.And(
                                pg.qubits[0] == swap.qubits[0],
                                pg.qubits[1] == swap.qubits[1],
                            ),
                            z3.And(
                                pg.qubits[0] == swap.qubits[1],
                                pg.qubits[1] == swap.qubits[0],
                            ),
                        )
                        or_lst.append(z3.And(to, eo, swap.on))
                    for qubit in pg.qubits:
                        for v in swap.qubits:
                            so = qubit.s == v.s  # Space overlap after wrapping
                            nc = z3.Not(z3.And(to, so, swap.on))  # No collision
                            and_lst.append(nc)

                eeo = z3.Or(
                    *or_lst
                )  # "Exists exact overlap" : there exists a swap that is on and exactly overlaps with the gate
                noo = z3.And(*and_lst)  # "No overlap at all"
                solver.add(
                    z3.Or(noo, eeo)
                )  # There is no overlap at all with an on swap gate, or there is an on swap gate that overlaps exactly with the gate.


def time(solver, phys_gates, depth):
    """
    All physical gates and routing swaps must be performed in the allotted circuit time. Note depth is a z3 var.
    """
    for pg in phys_gates.values():
        solver.add(0 <= pg.t, pg.t < depth)


def connectivity(solver, basis_graph, phys_gates):
    """
    Physical gates can only be performed on edges in a neighbourhood of 1 basis graph around the basis graph. Only edges with both of its vertices in a patch of 3 by 3 unit cells need to be taken into account.
    """

    def qubit_in_patch(q):
        return -1 <= q.x < 2 and -1 <= q.y < 2

    allowed_edges = basis_graph.get_open_patch(3, 3)
    allowed_edges = allowed_edges.translated(-1, -1)
    allowed_edges = list(allowed_edges.gates)
    allowed_edges = [
        edge
        for edge in allowed_edges
        if qubit_in_patch(edge.qubits[0]) and qubit_in_patch(edge.qubits[1])
    ]

    for pg in phys_gates.values():
        # if pg.n == 1:
        #     orlist = []
        #     for v in basis_graph.qubits:
        #         for i in range(-1, 2):
        #             for j in range(-1, 2):
        #                 orlist.append(pg.qubits[0] == v.translated(i, j))
        #     solver.add(z3.Or(*orlist))
        # elif pg.n == 2:
        if pg.n == 2:
            orlist = []
            for edge in allowed_edges:
                q0 = edge.qubits[0]
                q1 = edge.qubits[1]
                overlapsq = z3.Or(
                    z3.And(pg.qubits[0] == q0, pg.qubits[1] == q1),
                    z3.And(pg.qubits[0] == q1, pg.qubits[1] == q0),
                )
                orlist.append(overlapsq)
            solver.add(z3.Or(*orlist))


def consisteny(solver, phys_qubits, logic_qubits, phys_gates, depth):
    """
    Physical qubits should act on the physical qubits corresponding to the correct logical qubits.
    """
    r = len(phys_qubits)
    for lg, pg in phys_gates.items():
        for t in range(r):
            if lg.n == 1:
                lq = lg.qubits[0]
                phys_qubit_by_mapping = phys_qubits[t][lq]
                solver.add(
                    z3.Implies(
                        pg.t == t,
                        phys_qubit_by_mapping == pg.qubits[0],
                    )
                )
            elif lg.n == 2:
                lq0, lq1 = lg.qubits
                pq0_by_mapping = phys_qubits[t][lq0]
                pq1_by_mapping = phys_qubits[t][lq1]
                solver.add(
                    z3.Implies(
                        pg.t == t,
                        z3.And(
                            pq0_by_mapping == pg.qubits[0],
                            pq1_by_mapping == pg.qubits[1],
                        ),
                    )
                )


def cyclic(solver, phys_qubits):
    """
    If added, the final physical location of the logical qubits should be equal to the initial physical location of the logical qubits.
    """
    for lq, pq in phys_qubits[0].items():
        solver.add(phys_qubits[0][lq] == phys_qubits[-1][lq])


def gate_dependencies(solver, basis_circ, phys_gates):
    deps = basis_circ.get_gate_dependencies()
    for gate, gatep in deps:
        solver.add(phys_gates[gate].t < phys_gates[gatep].t)


def no_swap_swap_collisions(
    solver, rswaps
):  # Not superfluous. When there are as much logic qubits as there are physical qubits, swap-swap collisions are not possible because of injectivity of the logic-qubit-to-physical-qubit map `phys_qubits`. However, when there are much more physical qubits then there are logical qubits, two rSWAPs may be collide. Only needed I guess when the number of physical qubits is at least the number of logical qubits + 2.
    for T, rswapsT in enumerate(rswaps):
        for swap0, swap1 in it.combinations(rswapsT, 2):
            for i in range(2):
                for j in range(2):
                    solver.add(
                        z3.Implies(
                            swap0.qubits[i].s == swap1.qubits[j].s,
                            z3.Not(z3.And(swap0.on, swap1.on)),
                        )
                    )


def minimize_swaps(
    solver, rswaps, merge_swaps, phys_gates
):  # Only called when minimize_swaps==True
    if (
        merge_swaps == True
    ):  # If swaps can be merged, only penelize those swaps that cannot be merged.
        o_int_lst = []  # overlap int list
        for T, rswapsT in enumerate(rswaps):
            for swap in rswapsT:
                or_lst = []
                for pg in phys_gates.values():
                    to = pg.t == T  # Temporal overlap
                    for i in range(-1, 2):
                        for j in range(-1, 2):
                            eo = z3.Or(  # Exact spatial overlap
                                z3.And(
                                    pg.qubits[0].x + i == swap.qubits[0].x,
                                    pg.qubits[0].y + j == swap.qubits[0].y,
                                    pg.qubits[0].s == swap.qubits[0].s,
                                    pg.qubits[1].x + i == swap.qubits[1].x,
                                    pg.qubits[1].y + j == swap.qubits[1].y,
                                    pg.qubits[1].s == swap.qubits[1].s,
                                ),
                                z3.And(
                                    pg.qubits[0].x + i == swap.qubits[1].x,
                                    pg.qubits[0].y + j == swap.qubits[1].y,
                                    pg.qubits[0].s == swap.qubits[1].s,
                                    pg.qubits[1].x + i == swap.qubits[0].x,
                                    pg.qubits[1].y + j == swap.qubits[0].y,
                                    pg.qubits[1].s == swap.qubits[0].s,
                                ),
                            )
                            teo = z3.And(to, eo)  # The gate overlaps in space and time
                            or_lst.append(teo)
                egate = z3.Or(
                    *or_lst
                )  # There exists a gate that overlaps exactly in time and space (after translation) with the swap

                o_int_lst.append(
                    z3.If(z3.And(swap.on, z3.Not(egate)), 1, 0)
                )  # Add 1 to the list if the swap is on and does not overlap exactly with a gate.

        num_nak_swaps = z3.Int("num_nak_swaps")
        solver.add(num_nak_swaps == sum(o_int_lst))
        solver.minimize(num_nak_swaps)

    # Always minimize the total number of swap instructions.
    swap_vars = [z3.If(rswap.on, 1, 0) for rswapsT in rswaps for rswap in rswapsT]
    num_swap_vars = z3.Int("num_swap_vars")
    solver.add(num_swap_vars == sum(swap_vars))
    solver.minimize(num_swap_vars)
