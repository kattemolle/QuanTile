#!/usr/bin/env python3
"""
Main file containing definitions for QuanTile, the tilable quantum circuit compiler.
"""
from time import time
from copy import deepcopy
import z3
import networkx as nx
from datetime import datetime
import fcntl
import os
from . import constraints as cs

from .base import (
    Qubit,
    Gate,
    RSwap,
    Edge,
    BasisCirc,
    RoutedBasisCirc,
    BasisGraph,
    load_solution_database,
)


class Transpiler:
    """
    Route a `BasisCirc` to a device with the given `BasisGraph`.
    """

    def __init__(self, basis_circ: BasisCirc, basis_graph: BasisGraph):
        assert type
        self.basis_circ = basis_circ
        self.basis_graph = basis_graph
        self.name = "{} --> {}".format(self.basis_circ.name, self.basis_graph.name)
        self.solution = None  # Only populated after self.solve() is run.
        for gate in self.basis_circ.gates:
            assert (
                gate.name != "rSWAP"
            ), "Input circuit cannot contain routing swaps. Convert them to regular swaps, with e.g. name SWAP."

        # Transpiling options, to be set by directly setting the instance attribute. Below are the standard options, that you should override externally if you want to change them.
        self.cyclic = False
        self.fixed_depth = (
            False  # If set to an int, this int will be used as the fixed depth.
        )
        self.gate_dependencies = True
        self.merge_swaps = False
        self.slice_depth = None  # If set to an int d, the basis circuit is sliced into slices of depth d and the solver is run for these slices separately, after which the solutions are stitched together.
        self.init_map = None  # If set to a dict, take this dict as the initial map.
        self.final_map = None  # If set to a dict, take this dict as the final map.
        for gate in basis_circ.gates:
            assert (
                gate.name != "rSWAP"
            )  # Do not accept basis circuit that already have swaps in them FROM ROUTING.
            assert (
                gate.swap == False
            )  # Do not accept basis circuits that already have merged swaps in them
        self.sub_transpiler = False  # This is set to true manually if the transpiler is run as a sub transpiler originating from a non trivial slice_depth. Mainly for printing options.
        self.minimize_swaps = False  # After finding the lowest possible depth, minimize the number of swaps used at this depth.

    def __eq__(
        self, other
    ):  # Does not compare solutions, only names of input graphs and transpiling options.
        lst = [
            self.name == other.name,
            self.cyclic == other.cyclic,
            self.fixed_depth == other.fixed_depth,
            self.gate_dependencies == other.gate_dependencies,
            self.merge_swaps == other.merge_swaps,
            self.slice_depth == other.slice_depth,
            self.init_map == other.init_map,
            self.final_map == other.final_map,
            self.sub_transpiler == other.sub_transpiler,
            self.minimize_swaps == other.minimize_swaps,
        ]
        return all(lst)

    def solve(self):
        """
        Return a `BasisCirc` that implements `self.BasisCirc` while respecting the connectivity constraints of `self.BasisGraph`.
        """

        def declare_phys_seed_nums():
            """
            For the init map, we assume seed numbers are mapped to seed numbers. `phys_seeds[i]` gives the physical seed that seed i is mapped to.
            """
            # assert self.basis_circ.congruent_seeds()
            psns = {}
            for sq in self.basis_circ.get_seeds():
                ps = z3.Int("S_init_{}".format(sq.s))
                psns[sq.s] = ps
            return psns

        def declare_phys_qubits(depth, phys_seed_nums):
            """
            Declare and return z3 variables for all physical qubits that the logic_qubits are mapped to. The output `phys_qubits[T]` is a dict, where `phys_qubits[T][q]` is the phys `Qubit` that logical `Qubit` q is mapped to.
            """

            def physical_qubit_var(logic_qubit, T, phys_seed_nums):
                X = z3.Int(
                    "X_{}_{}_{}_{}".format(
                        logic_qubit.x, logic_qubit.y, logic_qubit.s, T
                    )
                )  # Note that the x,y,s,t values of the physical qubit are added to the X variable of the logical qubit to get a unique z3 variable name. That's the only reason for this name.
                Y = z3.Int(
                    "Y_{}_{}_{}_{}".format(
                        logic_qubit.x, logic_qubit.y, logic_qubit.s, T
                    )
                )
                if T == 0:
                    S = phys_seed_nums[logic_qubit.s]
                else:
                    S = z3.Int(
                        "S_{}_{}_{}_{}".format(
                            logic_qubit.x, logic_qubit.y, logic_qubit.s, T
                        )
                    )
                return Qubit(X, Y, S)

            phys_seed_nums = declare_phys_seed_nums()

            phys_qubits = []

            # A logic qubit can be mapped to any physical qubit. Note there can be rSWAPs in the last layer of the circ (in case of self.cyclic=True or if self.final_map!=None) so that phys_qubits needs depth+1 layers.
            r = depth
            if self.cyclic or self.final_map != None:
                r += 1
            for T in range(r):
                time_slice = {}
                for logic_qubit in self.basis_circ.qubits:
                    phys_qubit = physical_qubit_var(logic_qubit, T, phys_seed_nums)
                    time_slice[logic_qubit] = phys_qubit
                phys_qubits.append(time_slice)
            return phys_qubits

        def declare_phys_gates():
            """
            Declare and return dict with logic gates as keys and physical gates as values. The physical gates are populated by z3 vars. A physical gate has a name, qubits, and a time T at which it is placed.

            Note there are two types of physical qubits, the ones declared by declare_qubit_var and the ones declared by this function. In the function cs.consistency these two types are asserted to be consistent.
            """

            def phys_gate_var(logic_gate, i):
                phys_qubits = []
                for logic_qubit in logic_gate.qubits:
                    X = z3.Int(
                        "gate_X_x={}_y={}_s={}_t={}_index={}".format(
                            logic_qubit.x, logic_qubit.y, logic_qubit.s, logic_gate.t, i
                        )
                    )  # Note that the x,y,s,t,i values of the physical qubit are added to the X variable of the logical qubit to get a unique z3 variable name. This, and debugging, are the only reasons for this name.
                    Y = z3.Int(
                        "gate_Y_x={}_y={}_s={}_t={}__index={}".format(
                            logic_qubit.x, logic_qubit.y, logic_qubit.s, logic_gate.t, i
                        )
                    )
                    S = z3.Int(
                        "gate_S_x={}_y={}_s={}_t={}__index={}".format(
                            logic_qubit.x, logic_qubit.y, logic_qubit.s, logic_gate.t, i
                        )
                    )
                    phys_qubits.append(Qubit(X, Y, S))

                phys_qubits = tuple(phys_qubits)
                T = z3.Int(
                    "gate_T_x={}_y={}_s={}_t={}__index={}".format(
                        logic_gate.qubits[0].x,
                        logic_gate.qubits[0].y,
                        logic_gate.qubits[0].s,
                        logic_gate.t,
                        i,
                    )
                )  # This is just to get a unique variable name and for debugging. Note logic_gate.name and logic_gate.c are not z3 vars.
                return Gate(logic_gate.name, phys_qubits, T)

            phys_gates = {}
            for i, logic_gate in enumerate(self.basis_circ.gates):
                phys_gate = phys_gate_var(logic_gate, i)
                phys_gates[logic_gate] = phys_gate

            return phys_gates

        def declare_rswaps(depth):
            """
            Return the list rswaps, where rswaps[T] is the list of swaps at time T.
            """

            all_rswaps = []
            rswaps = []
            if self.cyclic or self.final_map != None:
                r = depth
            else:
                r = depth - 1
            for T in range(r):
                all_rswapsT = []  # All swaps at time T.
                rswapsT = []  # Only those swaps in the unit cell.
                for edge in self.basis_graph.edges:
                    on = z3.Bool(
                        "rSWAP_T={}_base_edge={}".format(T, edge.to_tuple())
                    )  # Only the edges in the basis graph get a new on variable.
                    rswap = RSwap(edge.qubits, T, on)
                    rswapsT.append(rswap)
                    for i in range(
                        -2, 3
                    ):  # A swap acting on the left side of a the basis circ at cell (2,0)  may swap a logical qubit from cell (1,0) to (2,0). This needs to be taken into account.
                        for j in range(-2, 3):
                            qubits = edge.translated(i, j).qubits
                            if (-1 <= qubits[0].x < 2 and -1 <= qubits[0].y < 2) or (
                                -1 <= qubits[1].x < 2 and -1 <= qubits[1].y < 2
                            ):
                                rswap = RSwap(qubits, T, on)
                                all_rswapsT.append(rswap)
                rswaps.append(rswapsT)
                all_rswaps.append(all_rswapsT)
            return all_rswaps, rswaps

        def phys_gates_model(phys_gates, model):
            """
            Return dict with keys the logical gates and values the physical gates that those are mapped to by routing.
            """
            phys_gates_model = {}
            for lg, pg in phys_gates.items():
                qubits_model = []
                for qubit in pg.qubits:
                    x_model = model[qubit.x]
                    y_model = model[qubit.y]
                    s_model = model[qubit.s]
                    assert (
                        type(x_model)
                        == type(y_model)
                        == type(s_model)
                        == z3.z3.IntNumRef
                    )
                    x_model = x_model.as_long()
                    y_model = y_model.as_long()
                    s_model = s_model.as_long()
                    qubit_model = Qubit(x_model, y_model, s_model)
                    qubits_model.append(qubit_model)
                qubits_model = tuple(qubits_model)
                t_model = model[pg.t].as_long()
                gate_model = Gate(pg.name, qubits_model, t_model)
                phys_gates_model[lg] = gate_model

            return phys_gates_model

        def swaps_model(swaps, model):
            swaps_model = []
            for T, swapsT in enumerate(swaps):
                for swap in swapsT:
                    if (
                        model[swap.on] == True
                    ):  # WARNING : this includes all on swaps from all cells to the model.
                        swap = Gate("rSWAP", swap.qubits, T)
                        swaps_model.append(swap)

            return swaps_model

        def map_model(phys_qubits, model):
            for T, apqT in enumerate(phys_qubits):
                for lq in apqT:
                    phys_qubits[T][lq].x = model[phys_qubits[T][lq].x].as_long()
                    phys_qubits[T][lq].y = model[phys_qubits[T][lq].y].as_long()
                    phys_qubits[T][lq].s = model[phys_qubits[T][lq].s].as_long()

            return phys_qubits

        def run_solver(depth, Smax):
            """
            Create and run the z3 solver at the given depth, return the z3 solver, the solution and the z3 vars used. If minimize_swaps==True, the number of swaps is reduced at the given depth.
            """
            if self.minimize_swaps == True:
                solver = z3.Optimize()
            else:
                solver = z3.Solver()

            # Declare z3 vars
            psns = declare_phys_seed_nums()
            pqs = declare_phys_qubits(depth, psns)
            pgs = declare_phys_gates()
            all_rswaps, rswaps = declare_rswaps(depth)

            lqs = self.basis_circ.qubits

            # Set constraints
            cs.mapping(
                solver,
                self.basis_graph,
                lqs,
                pqs,
                depth,
                Smax,
                self.init_map,
                self.final_map,
            )

            cs.time(solver, pgs, depth)
            cs.consisteny(solver, pqs, lqs, pgs, depth)
            cs.connectivity(solver, self.basis_graph, pgs)
            cs.swap_effect(
                solver, self.basis_graph, lqs, pqs, all_rswaps, self.init_map, Smax
            )
            if self.gate_dependencies == True:
                cs.gate_dependencies(solver, self.basis_circ, pgs)
            else:
                cs.no_gate_gate_collisions(solver, pgs)
            cs.no_swap_swap_collisions(solver, rswaps)
            cs.no_gate_swap_collisions(
                solver, pgs, all_rswaps, self.merge_swaps, rswaps
            )
            if self.cyclic:
                cs.cyclic(solver, pqs)
            if self.minimize_swaps:
                cs.minimize_swaps(
                    solver, rswaps, self.merge_swaps, pgs
                )  # TODO z3 Optimization module not efficient.

            sol = solver.check()

            return solver, sol, pqs, pgs, all_rswaps

        def post(solver, solved, start, end, depth, pqs, pgs, all_rswaps):
            """
            st
            """
            self.solution = {
                "solved": solved,
                "wall_clock": end - start,
                "depth": depth,
                "map_model": None,
                "routed_basis_circ": None,
            }

            if solved == True:
                mod = solver.model()
                apgm = phys_gates_model(pgs, mod).values()
                sm = swaps_model(all_rswaps, mod)

                # If self.merge_swaps == True, merge any swaps that overlap EXACTLY with the gate into that gate by setting gate.swap=True and remove equivalent swaps in other cells. Swaps are to be performed after the gate.
                if self.merge_swaps == True:
                    _sm = deepcopy(sm)
                    for gate in apgm:
                        for swap in _sm:
                            so = (
                                gate.qubits == swap.qubits
                                or gate.qubits == swap.qubits[::-1]
                            )
                            to = gate.t == swap.t
                            if so and to:
                                gate.swap = True
                                for i in range(-2, 3):
                                    for j in range(-2, 3):
                                        tr = swap.translated(i, j)
                                        if tr in sm:
                                            sm.remove(tr)

                # Remove any swaps sill left that are not in the basis graph.
                sm = [s for s in sm if Edge(s.qubits) in self.basis_graph.edges]
                # print(self.basis_graph.edges)

                apgm = list(apgm) + sm
                routed_name = "{}_on_{}".format(
                    self.basis_circ.name, self.basis_graph.name
                )
                apgm.sort(key=lambda gate: gate.t)

                mm = map_model(pqs, mod)
                rbc = RoutedBasisCirc(
                    apgm, mm, name=routed_name
                )  # To route a circ, it is assumed it has congruent seeds. The  routed circuit may not have congruent seeds, so it may not be possible to route a routed circuit again, but this is not a use case anyway.

                if self.fixed_depth == None:
                    for layer in rbc.get_layers():
                        assert (
                            len(layer) != 0
                        ), "In an optimal-depth circuit there can be no empty layers."

                self.solution["routed_basis_circ"] = rbc
                self.solution["map_model"] = mm

        # End of defs in solve().
        assert (
            self.basis_graph.get_smax() >= self.basis_circ.get_smax()
        ), "The number of seeds in the basis graph must at least be the number of seed qubits in the basis circuit. Take a bigger patch, reseed it, and set it as the new input basis graph."

        print()
        if self.sub_transpiler == False:
            print("Solving")
            print(vars(self), flush=True)
        else:
            print("    Sub-solving")
            print("   ", vars(self), flush=True)

        if self.slice_depth != None:
            return self.sliced_solve()

        solved = False
        self.solution = {"time_stamp": datetime.now()}
        if type(self.fixed_depth) == int:
            self.init_depth = self.fixed_depth
        else:
            if self.gate_dependencies == True:
                self.init_depth = self.basis_circ.get_critical_path_length()
            else:
                _bc = self.basis_circ.get_open_patch(3, 3)
                self.init_depth = _bc.get_max_qubit_activity()

        depth = deepcopy(self.init_depth)
        Smax = self.basis_graph.get_smax()
        try:
            db = load_solution_database()
            for group in db:
                if self in db[group]:
                    print(
                        "WARNING: transpiler already in database at group '{}', entry {}. Transpiling nevertheless (not adding to database by default.)".format(
                            group, db[group].index(self)
                        )
                    )
        except:
            print("Failed loading database with previous solutions.")

        start = time()
        while not solved:
            if self.sub_transpiler == False:
                print("trying depth =", depth, "...")
            else:
                print("    trying depth =", depth, "...")

            (solver, sol, pqs, pgs, all_rswaps) = run_solver(depth, Smax)

            if type(self.fixed_depth) == int:
                if sol == z3.sat:
                    if self.sub_transpiler == False:
                        print("solved")
                    else:
                        print("    solved")
                    solved = True
                else:
                    print("  unsat,  but stopping search")
                    break
            else:
                if sol == z3.sat:
                    solved = True
                    if self.sub_transpiler == False:
                        print("solved")
                    else:
                        print("    solved")
                else:
                    depth += 1

        end = time()

        post(
            solver, solved, start, end, depth, pqs, pgs, all_rswaps
        )  # Postprocess and store the solution in self.solution

        return self.solution

    def sliced_solve(self):
        """
        Schedule all gates as early as possible and slice up the circuit in slices with depth d. Solve the routing problem optimally for these subcircuits.
        """
        assert self.slice_depth != None
        assert type(self.slice_depth) == int
        assert self.slice_depth > 0

        bc = self.basis_circ.rescheduled()
        layers = bc.get_layers()
        bcs = []
        for i in range(0, len(layers), self.slice_depth):
            bc = layers[i : i + self.slice_depth]
            bc = [gate for layer in bc for gate in layer]
            bc = BasisCirc(bc, qubits=self.basis_circ.qubits)
            bcs.append(bc)

        ts = []
        for i, bc in enumerate(bcs):
            t = Transpiler(bc, self.basis_graph)
            t.name += "slice_{}".format(i)
            t.cyclic = False
            t.gate_dependencies = self.gate_dependencies
            t.merge_swaps = self.merge_swaps
            t.sub_transpiler = True
            t.slice_depth = None
            ts.append(t)

        assert (
            len(ts) > 1
        ), "No slicing needed with current parameters, just use solve()"

        ts[0].solve()
        i = 1
        while i < len(ts) - 1:
            final_map = ts[i - 1].solution["map_model"][-1]
            ts[i].init_map = final_map
            ts[i].solve()
            i += 1

        final_map = ts[i - 1].solution["map_model"][-1]
        ts[i].init_map = final_map
        if self.cyclic == True:
            ts[i].final_map = ts[0].solution["map_model"][0]
        ts[i].solve()

        # Stitch together solutions
        self.solved = True

        # Shift all times of gates by the correct amount and add up wall clock times
        wall_clock = 0
        i = 0
        tmax = 0
        rbcs = []
        while i < len(ts):
            t = ts[i]
            wall_clock += t.solution["wall_clock"]
            rbc = t.solution["routed_basis_circ"]

            if i != 0:
                for gate in rbc.gates:
                    gate.t += tmax + 1

            rbcs.append(rbc)
            tmax = rbc.get_tmax()
            i += 1

        # Add boundary swaps to gates list
        gates = [gate for rbc in rbcs for gate in rbc.gates]

        # Condense the routed basis circuit
        rbc = BasisCirc(gates).rescheduled()
        gates = rbc.gates

        # Add boundary swaps to gates list
        egates = []
        for gate in gates:
            if gate.name == "rSWAP":
                for i in range(-2, 3):
                    for j in range(-2, 3):
                        egates.append(gate.translated(i, j))
            elif gate.swap == True:
                for i in range(-2, 3):
                    for j in range(-2, 3):
                        if not (i == 0 and j == 0):
                            rSWAP = Gate("rSWAP", gate.qubits, t=gate.t)
                            egates.append(rSWAP.translated(i, j))
                        else:
                            egates.append(gate)
            else:
                egates.append(gate)

        erbc = BasisCirc(
            egates,
            collision_checking=False,
            unique_seeds_checking=False,
            wrapped_circ_checking=False,
        )  # Enlarged routed basis circuit. It is completed with the boundary rSWAPs. TODO : can this be replaced by get_completed_patch()?

        # Reconstruct map model from enlarged routed basis circ.
        layers = erbc.get_layers()
        init_map = ts[0].solution["map_model"][0]
        map_model = [init_map]
        for layer in layers:
            map_layer = deepcopy(map_model[-1])
            for gate in layer:
                if gate.name == "rSWAP" or gate.swap == True:
                    pqs = gate.qubits
                    inv_map = {v: k for k, v in map_layer.items()}
                    if pqs[0] in inv_map.keys():
                        lq0 = inv_map[pqs[0]]
                        map_layer[lq0] = pqs[1]
                    if pqs[1] in inv_map.keys():
                        lq1 = inv_map[pqs[1]]
                        map_layer[lq1] = pqs[0]

            map_model.append(map_layer)

        assert map_model[-1] == ts[-1].solution["map_model"][-1]

        rbc = RoutedBasisCirc(rbc.gates, map_model)

        self.solution = {
            "solved": True,
            "wall_clock": wall_clock,
            "depth": rbc.get_tmax() + 1,
            "map_model": map_model,
            "routed_basis_circ": rbc,  # Note this is the routed basis circuit without the boundary swaps
        }

        print("solved")

        return self.solution

    def plot_solution(self, folder="plots"):
        """
        Put a sequence of pdfs in folder `folder` that depict the layers of the circuit. A 3x3 patch of the device connectivity graph is depicted with grey edges. Blue edges depict 2-qubit gates from the circuit. Blue edges depict the swaps. The solution basis circuit is bold.
        """
        assert self.solution != None, "First run self.solve()"

        pyplotc = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]
        otherc = ["#2f17a9", "#99f14a"]
        colors = pyplotc + otherc
        colors = [color + "99" for color in colors]  # Set alpha

        ldg = self.basis_graph.get_plot_graph()

        c = self.solution["routed_basis_circ"]
        c = c.get_open_patch(2, 2, -1, -1)
        layers = c.get_layers()

        for l, layer in enumerate(layers):
            g = deepcopy(ldg)
            g = nx.MultiGraph(
                g
            )  # Allow multi-edges for also plotting input circuit connectivity graph.
            for gate in layer:
                if gate.n == 2:  # Single-qubit gates are not plotted.
                    edge = (
                        gate.qubits[0].to_tuple(),
                        gate.qubits[1].to_tuple(),
                    )
                    if gate.name == "rSWAP":
                        g[edge[0]][edge[1]][0]["color"] = colors[1]
                    elif gate.swap == False:
                        g[edge[0]][edge[1]][0]["color"] = colors[0]
                    else:
                        g[edge[0]][edge[1]][0]["color"] = colors[5]

            import warnings

            warnings.filterwarnings(
                "ignore", category=RuntimeWarning, module="pygraphviz"
            )

            g = nx.nx_agraph.to_agraph(g)
            g.layout(prog="neato", args="-Goutputorder=edgesfirst -Gsplines=true")

            # Add connectivity graph of input circuit
            con_graph = self.basis_circ.to_basis_graph()
            # con_graph = con_graph.get_patch(3, 3, 1)[0]
            # con_graph = con_graph.translated(-1, -1)

            mm = self.solution["map_model"]
            for edge in con_graph.gates:
                v0, v1 = edge.qubits
                pq0 = mm[l][v0]
                pq1 = mm[l][v1]
                g.add_edge(pq0, pq1, color=colors[4], key=1)

            if not os.path.exists(folder):
                os.mkdir(folder)

            g.draw("{}/{}_{}.png".format(folder, self.name, str(l).zfill(4)))

    def add_to_database(self, group="unsrt", fname="solutions.pkl"):
        """
        Add transpiler (with possibly solution) to the database under the group name.
        """
        import pickle

        if os.path.isfile(fname):
            with open(fname, "rb") as f:
                db = pickle.load(f)
        else:
            db = {}

        if group in db and self in db[group]:
            print(
                "WARNING: transpiler already in database. Adding possibly duplicate entry."
            )

        if group not in db:
            db[group] = []

        db[group].append(self)

        with open(fname, "wb") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            pickle.dump(db, f)
            fcntl.flock(f, fcntl.LOCK_UN)

    def append_to_qasm_database(self, fname="solutions_qasm.txt"):
        assert self.solution != None
        qasm = self.solution["routed_basis_circ"].to_qasm()

        with open(fname, "a") as f:
            f.write(str(vars(self)) + "\n")
            f.write("// begin of qasm file for the routed basis circuit\n")
            f.write(qasm + "\n")
            f.write("// end of qasm file for the routed basis circuit\n")
