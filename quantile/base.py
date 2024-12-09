#!/usr/bin/env python3

import numpy as np
import json
from copy import copy, deepcopy
import itertools as it
import z3
import ast
import networkx as nx
import os
from importlib.resources import files


class Qubit:
    """
    Class to define a qubit.

    Parameters:
    -----------
      x : int or
        x coordinate of cell the qubit is in
      y : int
        y coordinate
      s : int
        seed number
    """

    def __init__(self, x, y, s):
        assert (type(x) == type(y) == type(s) == int) or (
            type(x) == type(y) == type(s) == z3.z3.ArithRef
        )
        self.x = x
        self.y = y
        self.s = s
        self.flag = False

    def __eq__(self, other):
        if type(other) == Qubit:
            if type(self.x) == type(self.y) == type(self.s) == int:
                return self.to_tuple() == other.to_tuple()
            elif type(self.x) == type(self.y) == type(self.s) == z3.z3.ArithRef:
                return z3.And(self.x == other.x, self.y == other.y, self.s == other.s)
            else:
                return TypeError
        else:
            return False

    def __ne__(self, other):
        if type(self.x) == type(self.y) == type(self.s) == int:
            return not self.__eq__(other)
        elif type(self.x) == type(self.y) == type(self.s) == z3.z3.ArithRef:
            return z3.Not(self.__eq__(other))
        else:
            return TypeError

    def __lt__(self, other):
        return self.to_tuple() < other.to_tuple()

    def __str__(self):
        return str(self.to_tuple())

    def __repr__(self):
        return "<Qubit " + self.__str__() + ">"

    def __hash__(self):
        return hash(self.to_tuple())

    def translated(self, dx, dy):
        """
        Return qubit that is translated dx cells horizontally and dy cells vertically.
        """
        qubit = deepcopy(self)
        qubit.x += dx
        qubit.y += dy
        return qubit

    def to_tuple(self):
        return (self.x, self.y, self.s, self.flag)


class Gate:
    """
    Class to define a gate.
    """

    def __init__(self, name, qubits, t=None):
        if type(t) == int:
            assert t >= 0
        self.name = name  # If name is in the list  ["p", "x", "y", "z", "h", "s", "sdg", "t", "tdg", "sx", "rx", "rz", "cx", "cy", "cz", "sp", "crx", "cry", "crz", "ch", "cu", "swap", "U", "gphase"]  it will be interpreted as a pre-defined QASM gate upon export to QASM. It will overwrite any custom definition of those gates in self.qasm_unitary.
        self.n = len(qubits)
        assert self.n in [1, 2]
        if self.n == 2:
            for q, qp in it.combinations(qubits, 2):
                if type(q.s) == int:  # Only checks the type of q.s...
                    assert not (q == qp), "A gate cannot act on the same qubit twice"
        self.qubits = qubits
        self.t = t
        self.swap = False  # If True, the gate is to be interpreted as that gate, FOLLOWED by a swap.
        if self.n == 1:
            self.qasm_unitary = "id a;"  # Alter to make custom literal QASM definition of the unitary. If `gate.name` is in the OpenQASM 3 standard gate library or U or gphase, (https://openqasm.com/language/standard_library.html#standard-library), this definition (inculding the current identity definition) is overruled upon export to OpenQASM. Currently only one and two-qubit gates are supported. Currently, you always have to define all gates with a custom qasm unitary if merge_swaps==True in the transpiler and you want to export the result to qams with correct unitary definitions.
        elif self.n == 2:
            self.qasm_unitary = "id a;\n  id b;"  # Same as for single qubit gates.
        else:
            raise ValueError

    def __eq__(self, other):
        if type(other) == Gate:
            return (
                self.name == other.name
                and self.qubits == other.qubits
                and self.t == other.t
            )
        else:
            return False
        # Warning : self.pre_wrapped, if set at all, is not taken into account.

    def __lt__(self, other):
        return hash(self) < hash(other)

    def __str__(self):
        return "{} {} {} {}".format(self.name, self.qubits, self.t, self.swap)
        # Warning : self.pre_wrapped, if set at all, is not taken into account.

    def __repr__(self):
        return "<Gate " + self.__str__() + ">"
        # Warning : self.pre_wrapped, if set at all, is not taken into account.

    def __hash__(self):
        return hash((self.name, self.qubits, self.t))
        # Warning : self.pre_wrapped, if set at all, is not taken into account.

    def translated(self, dx, dy, dt=0):
        """
        Return translated gate dx unit cells horizontally, dy unit cells vertically and dt in time.
        """
        gate = deepcopy(self)
        gate.qubits = tuple(qubit.translated(dx, dy) for qubit in self.qubits)
        if gate.t != None:
            gate.t += dt
        return gate

    def acts_on_logical_qubit(self, phys_qubits):
        """
        Return `True` if the gate acts on a logical qubit. `phys_qubits[t]` is a dict with logical qubits at time t as keys.
        """
        for pq in self.qubits:
            if pq in phys_qubits[self.t].values():
                return True
        return False


class Circ:
    """
    Circuit. Basically as list of `Gate`s, where every qubit in every Gate is a `Qubit`, with x,y coordinates and a seed.
    Users usually interact with `BasisCirc`, not directly with `Circ`.
    """

    def __init__(self, gates, qubits=None, name=None, collision_checking=True):
        assert type(gates) == list
        self.gates = gates
        self.name = name
        if collision_checking:
            assert self.no_collisions()
        self.qubits = qubits
        if qubits is None:
            self.qubits = []
        self.qubits = self.get_active_qubits().union(set(self.qubits))
        self.n = len(self.qubits)
        self.collision_checking = collision_checking

    def __str__(self):
        return str(self.gates)

    def get_tmax(self):
        """
        Latest time at which a gate acts.
        """
        tmax = 0
        for gate in self.gates:
            if gate.t != None and gate.t > tmax:
                tmax = gate.t

        return tmax

    def get_layers(self):
        assert all(
            gate.t != None for gate in self.gates
        ), "first assign times to the gates to get the layers"
        depth = self.get_tmax() + 1
        if type(self.gates[0]) == Edge:
            ls = [self.gates]
        else:
            ls = [[] for _ in range(depth)]
            for gate in self.gates:
                ls[gate.t].append(gate)
        ls = deepcopy(ls)
        return ls

    def get_all_active_qubits(self):
        """
        Return a LIST of all `Qubit`s that the circuit acts on. May contain duplicates.
        """
        all_qubits = []
        for gate in self.gates:
            for qubit in gate.qubits:
                all_qubits.append(qubit)

        all_qubits = deepcopy(all_qubits)
        return all_qubits

    def get_active_qubits(self):
        """
        Return the SET of all `Qubit`s (no duplicates).
        """
        return set(self.get_all_active_qubits())

    def get_seed_nums(self):
        """
        Return a list of seed numbers of the `Qubit`s that the circuit acts on (no duplicates).
        """
        all_seed_nums = [qubit.s for qubit in self.qubits]
        return list(set(all_seed_nums))

    def get_smax(self):
        """
        Return the largest seed number (int) in any of the `Qubit`s that the `Circ` acts on.
        """
        return max(self.get_seed_nums())

    def get_max_qubit_activity(self):
        activity = {}
        for gate in self.gates:
            for qubit in gate.qubits:
                if qubit in activity:
                    activity[qubit] += 1
                else:
                    activity[qubit] = 1
        vals = activity.values()
        return max(vals)

    def translated(self, dx, dy, dt=0):
        """
        Return circuit copy, translated by dx cells horizontally, dy cells vertically, and dt in time.
        """
        circ = deepcopy(self)
        circ.gates = list(gate.translated(dx, dy, dt) for gate in self.gates)
        return circ

    def no_collisions(self):
        for i, layer in enumerate(self.get_layers()):
            qubits = []
            for gate in layer:
                for qubit in gate.qubits:
                    qubits.append((qubit.x, qubit.y, qubit.s))
            if len(qubits) != len(set(qubits)):
                import pickle

                with open("dump_no_col.pkl", "wb") as f:
                    d = (self, i, layer, qubits, len(qubits), len(set(qubits)))
                    pickle.dump(d, f)
                print("Collision encountered, dump completed.")
                return False
        return True

    def repeated(self, nc):
        """
        Repeat a copy of self repeated temporally nc times. Does not also reschedule.
        """
        if nc == 1:
            return self
        else:
            gates = self.gates
            depth = self.get_tmax() + 1
            cycles = [deepcopy(gates) for _ in range(nc)]
            for cy, _ in enumerate(cycles):
                for g, _ in enumerate(cycles[cy]):
                    cycles[cy][g].t = cycles[cy][g].t + cy * depth
            gates = [gate for cycle in cycles for gate in cycle]
            if type(self) == Circ:
                c = Circ(gates)
            if type(self) == BasisCirc:
                c = BasisCirc(
                    gates,
                    self.name,
                    self.collision_checking,
                    self.unique_seeds_checking,
                    self.congruent_seeds_checking,
                    self.wrapped_circ_checking,
                    self.qubits,
                )
            return c

    def to_dag(self):
        """
        Convert the `Circ` circ to the directed acyclic graph (DAG) representation, which is outputted as a nx.DiGraph, with `Gate` objects as nodes and `Qubit` objects as edge labels. Removes idling qubits.
        """
        g = nx.MultiDiGraph()
        qubits = self.get_active_qubits()

        for qubit in qubits:
            g.add_edge((qubit, "init"), qubit)

        for layer in self.get_layers():
            for gate in layer:
                g.add_node(gate)
                for qubit in gate.qubits:
                    assert qubit in g.nodes(), f"{qubit},{g.nodes}"
                    pre = list(g.predecessors(qubit))
                    npre = len(pre)
                    assert npre in [0, 1], (qubit, pre, len(pre))
                    if npre == 1:
                        pre = pre[0]
                        assert (pre, qubit, 1) not in g.edges(keys=True)
                        nx.set_edge_attributes(g, {(pre, qubit, 0): {"qubit": qubit}})
                    g = nx.contracted_nodes(g, gate, qubit)
                    g.add_edge(gate, qubit, qubit=qubit)

        mapping = {n: (n, "final") for n in g if type(n) == Qubit}
        gp = deepcopy(g)  # Seems to circumvent a bug in networkx
        g = nx.relabel_nodes(gp, mapping)

        assert nx.is_directed_acyclic_graph(
            g
        ), "Graph representation of the circuit is not a directed acyclic graph."

        return g

    @classmethod
    def from_dag(cls, dag, qubits=None, name=None):
        """
        Construct self from a DAG circuit `dag` (nx.MultiDiGraph). Ignores time coordinates of gates. Idling qubits can be added with the qubits argument.
        """

        def collision(gate1, gate2):
            if gate1.t == gate2.t:
                for q in gate1.qubits:
                    for p in gate2.qubits:
                        if q == p:
                            return True
            return False

        def remove_gate(dag, gate):
            _dag = deepcopy(dag)
            ies = sorted(
                list(_dag.in_edges(gate, data=True, keys=True)),
                key=lambda e: e[3]["qubit"].to_tuple(),
            )  # In edges, sorted by qubit label
            oes = sorted(
                list(_dag.out_edges(gate, data=True, keys=True)),
                key=lambda e: e[3]["qubit"].to_tuple(),
            )  # Out edges, sorted by qubit label
            assert gate in _dag
            edges = list(_dag.edges(keys=True, data=True))
            _dag = nx.MultiDiGraph(
                [
                    edge for edge in edges if edge[0] != gate and edge[1] != gate
                ]  # Networkx implementation for node removal has a bug...
            )
            q0 = ies[0][3]["qubit"]
            assert q0 == oes[0][3]["qubit"]
            q0p = deepcopy(q0)
            _dag.add_edge(ies[0][0], oes[0][1], qubit=q0p)
            if gate.n == 2:
                q1 = ies[1][3]["qubit"]
                assert q1 == oes[1][3]["qubit"]
                q1p = deepcopy(q1)
                _dag.add_edge(ies[1][0], oes[1][1], qubit=q1p)
            return deepcopy(_dag)

        assert type(dag) == nx.MultiDiGraph
        _dag = deepcopy(dag)
        gates = []
        clean = False
        t = 1
        while not clean:
            nodes = iter(list(_dag.nodes))
            while True:
                try:
                    node = next(nodes)
                except StopIteration:
                    clean = True
                    break
                pre = list(_dag.predecessors(node))
                if (
                    type(node) == Gate
                    and len(pre) != 0
                    and all(type(m) == tuple and m[1] == "init" for m in pre)
                ):
                    gate = deepcopy(node)
                    gate.t = t
                    col = False
                    while not col and gate.t >= 0:
                        gate.t -= 1
                        col = any(collision(gate, sgate) for sgate in gates)
                    gate.t += 1
                    gates.append(gate)
                    assert node in _dag
                    _dag = remove_gate(_dag, node)
                    t += 1
                    break

        return cls(
            gates,
            qubits=qubits,
            name=name,
        )

    def to_wrapped_circ(self):
        raise NotImplementedError("Wrapping only supported for `BasisGraph`s")

    def to_qasm(self):
        """
        Return an OpenQASM string. Be sure that all gate names comply with OpenQASM rules for gate
        """
        circ = deepcopy(self)
        circ.gates.sort(key=lambda x: x.t)
        qasm = ""

        cgns = []  # custom declared gate names
        qgns = [
            "p",
            "x",
            "y",
            "z",
            "h",
            "s",
            "sdg",
            "t",
            "tdg",
            "sx",
            "rx",
            "rz",
            "cx",
            "cy",
            "cz",
            "sp",
            "crx",
            "cry",
            "crz",
            "ch",
            "cu",
            "swap",
            "U",
            "gphase",
        ]  # predefined 1- and 2-q OpenQASM gate names.

        def replace_min(q):  # "-" is not allowed in qasm... so replace with min
            q = deepcopy(q)
            if q.x < 0:
                q.x = "min" + str(abs(q.x))
            if q.y < 0:
                q.y = "min" + str(abs(q.y))
            return q

        def gate_replace_min(g):
            if g.n == 1:
                q = g.qubits[0]
                q = replace_min(q)
                g.qubits = (q,)
            if g.n == 2:
                q0, q1 = g.qubits
                q0 = replace_min(q0)
                q1 = replace_min(q1)
                g.qubits = (q0, q1)

            return g

        qasm += "OPENQASM 3;\n"
        qasm += 'include "stdgates.inc";\n'
        qasm += "\n"

        # Declare qubits
        qasm += "// declare qubits\n"
        for q in circ.qubits:
            q = replace_min(q)
            qasm += f"qubit q_{q.x}_{q.y}_{q.s};\n"
        qasm += "\n"

        # Declare gates
        qasm += "// declare gates\n"

        for gate in circ.gates:
            gate_replace_min(gate)

        for gate in circ.gates:
            if gate.swap == True:
                gate.name = "swap_" + gate.name
                gate.qasm_unitary += "\n  swap a,b;"

            split = gate.name.split("(")
            gname = split[0]  # Gate name before parameter

            if len(split) > 1:
                par = "(la)"
            else:
                par = ""

            if gname not in qgns and gname not in cgns and gname != "rSWAP":
                cgns.append(gname)
                if gate.n == 1:
                    qasm += f"gate {gname}{par} a\n{{"
                elif gate.n == 2:
                    qasm += f"gate {gname}{par} a, b\n{{"
                else:
                    raise ValueError
                qasm += f"\n  {gate.qasm_unitary}\n}}\n"
        qasm += "\n"

        # Write gates:
        qasm += "// circuit\n"
        for gate in circ.gates:
            if gate.n == 1:
                q = gate.qubits[0]
                qasm += f"{gate.name} q_{q.x}_{q.y}_{q.s};\n"
            elif gate.n == 2:
                if gate.name == "rSWAP":
                    gname = "swap"  # This is the only point where rSWAPs are cast to normal swaps.
                else:
                    gname = gate.name  # includes the parameter as a float
                q0, q1 = gate.qubits
                qasm += f"{gname} q_{q0.x}_{q0.y}_{q0.s}, q_{q1.x}_{q1.y}_{q1.s};\n"
            else:
                raise ValueError

        return qasm

    def export_qasm(self):
        """
        Export to qasm to self.name.qasm.
        """
        if self.name == None:
            name = "circ"
        else:
            name = self.name
        with open(f"{name}.qasm", "w") as f:
            f.write(self.to_qasm())


class BasisCirc(Circ):
    """
    Basis circuit, the 'tile' of a tilable circuit. The set of seed numbers of the qubits that the `Circ` acts on must be of the form {0,1,2,3,...nS-1}.
    """

    def __init__(
        self,
        gates,
        name=None,
        collision_checking=True,
        unique_seeds_checking=True,
        congruent_seeds_checking=True,
        wrapped_circ_checking=True,
        qubits=None,
    ):
        if qubits is None:
            self.qubits = []

        if all(gate.t == None for gate in gates):
            t = 0
            gates = deepcopy(gates)
            for gate in gates:
                gate.t = t
                t += 1
        else:
            for gate in gates:
                assert type(gate.t) == int

        super().__init__(
            deepcopy(gates),
            qubits=qubits,
            name=name,
            collision_checking=collision_checking,
        )
        self.unique_seeds_checking = unique_seeds_checking
        self.congruent_seeds_checking = congruent_seeds_checking
        self.wrapped_circ_checking = wrapped_circ_checking
        self.phys_qubits = None

        if unique_seeds_checking:
            assert self.unique_seeds()
        if congruent_seeds_checking:
            assert self.congruent_seeds(), sorted(self.get_seed_nums())
        if wrapped_circ_checking == True:
            self.to_wrapped_circ()  # This checks that the wrapped circ is a valid Circ.

    def unique_seeds(self):
        """
        Return `True` if every seed occurs in every layer of `self.circ` at most once.
        """
        for layer in self.get_layers():
            seeds = []
            for gate in layer:
                for qubit in gate.qubits:
                    seeds.append(qubit.s)
            if len(seeds) != len(set(seeds)):
                return False
        return True

    def congruent_seeds(self):
        """
        Return `True` if the set of seed numbers is of the form {0,1,2,3,...nS-1}.
        """
        smax = self.get_smax()
        return sorted(self.get_seed_nums()) == list(range(smax + 1))

    def get_seeds(self):
        """
        Return the SET of qubits that are in the unit cell (no duplicates).
        """
        ucq = {q for q in self.qubits if (q.x == 0 and q.y == 0)}
        return deepcopy(ucq)

    def get_gates(self):
        return self.gates

    def get_open_patch(self, n, m, ln=0, lm=0, nc=1):
        """
        Return circuit on a patch of n (horizontal) by m (vertical) by nc (time direction) circuit cells.
        `lm` and `ln` set the lower bound of the coordinate. If set, the basis graph's indices range over range(ln,n) and range(lm,m).
        """
        gates = []
        for dx in range(ln, n):
            for dy in range(lm, m):
                circ = self.translated(dx, dy)
                gates += circ.gates

        if type(self) == BasisGraph:
            assert len(gates) == len(
                set(gates)
            ), "Basis Graph has redundant edges. Please remove any redundant edge."

            p = Circ(gates, collision_checking=False)
        else:
            p = Circ(gates, collision_checking=True).repeated(nc)

        return p

    def get_completed_patch(self, n, m, phys_qubits=None, nc=1):
        """
        Return circuit on a patch of n (horizontal) by m (vertical) by nc (time direction) circuit cells completed with the appropriate extra rSWAPs coming from other tiles.

        For the format of phys_qubits, see the format of map_model in transpiler.solve()

        For a routed basis circ to be equivalent to the basis circ, boundary rSWAPs need to be taken into account. So even when we route a single basis circuit, that routed circuit becomes equivalent to the basis circuit only after calling get_completed_patch(1,1,phys_qubits)

        """
        assert type(self) != BasisGraph
        assert phys_qubits != None
        if nc > 1:
            assert (
                phys_qubits[0] == phys_qubits[-1]
            ), "Repeating this circuit is not compatible with the logical-to-physical qubit map."

        # Enlarge set of phys qubits
        epq = []
        assert phys_qubits != None
        for time_slice in phys_qubits:
            ts = {}
            for lq, pq in time_slice.items():
                for i in range(n):
                    for j in range(m):
                        ts[lq.translated(i, j)] = pq.translated(i, j)
            epq.append(ts)

        # Construct the boundary
        lp = self.get_open_patch(n + 3, m + 3, -2, -2).gates
        sp = self.get_open_patch(n, m).gates
        bdry = list(set(lp) - set(sp))
        _bdry = []
        for g in bdry:
            if g.name == "rSWAP" and g.acts_on_logical_qubit(epq):
                _bdry.append(g)
            elif g.swap == True and g.acts_on_logical_qubit(epq):
                rSWAP = Gate("rSWAP", g.qubits, t=g.t)
                _bdry.append(rSWAP)
        bdry = _bdry

        # Merge boundary and the patch
        p = sp + bdry
        p = Circ(p)
        p = p.repeated(nc)
        # p = p.rescheduled()

        return p

    get_patch = get_open_patch

    def to_basis_graph(self):
        """
        Return connectivity graph that the circuit uses.
        """
        assert type(self) == BasisCirc
        edges = []
        for gate in self.gates:
            edge = Edge(gate.qubits)
            edges.append(edge)

        edges = deepcopy(edges)
        return BasisGraph(edges, name=self.name)

    def resized(self, n, m):
        """
        Get a patch of the lattice circuit that is induced by the basis circuit (of n by m basis cirucits) and reseed it. It gives a new basis circuit that is equivalent to a circuit of n by m old basis circuits. Also returns the mapping applied to all qubits in all gates. The name of the new basis circ contains the new size n,m.
        """
        patch_gates = self.get_open_patch(n, m).gates
        nc = 1  # Number of cycles
        assert n > 0 and m > 0 and nc >= 1
        seed_qubits = {
            q
            for gate in patch_gates
            for q in gate.qubits
            if 0 <= q.x < n and 0 <= q.y < m
        }

        mapping = {}
        for nsn, q in enumerate(seed_qubits):
            mapping[q] = Qubit(0, 0, nsn)

        nonseed_qubits = {
            q
            for gate in patch_gates
            for q in gate.qubits
            if not (0 <= q.x < n and 0 <= q.y < m)
        }  # This uses that all seed qubits of the new basis circuit have an edge connected to them.
        for nsq in nonseed_qubits:
            nsn = mapping[Qubit(nsq.x % n, nsq.y % m, nsq.s)].s
            mapping[nsq] = Qubit(
                int(np.floor(nsq.x / n)), int(np.floor(nsq.y / m)), nsn
            )

        gates = []
        for gate in patch_gates:
            nq0 = mapping[gate.qubits[0]]
            nqs = (nq0,)
            if gate.n == 2:
                nq1 = mapping[gate.qubits[1]]
                nqs = (nq0, nq1)
            gate.qubits = nqs
            gates.append(gate)

        old_gates = [
            (gate.qubits[0], gate.qubits[1]) for gate in patch_gates if gate.n == 2
        ]
        new_gates = [(gate.qubits[0], gate.qubits[1]) for gate in gates if gate.n == 2]
        patch_graph = nx.Graph(old_gates)
        new_patch_graph = nx.Graph(new_gates)
        assert nx.is_isomorphic(patch_graph, new_patch_graph)

        if self.name != None:
            name = self.name + f"({n},{m})"
        else:
            self.name = f"({n},{m})"

        if type(self) == BasisCirc:
            return BasisCirc(gates, name=name), mapping
        elif type(self) == BasisGraph:
            return BasisGraph(gates, name=name), mapping

    def to_wrapped_circ(self):
        """
        Wrap the qubits of a deepcopy around the edges of the unit cell and return as a Circ. Wrapping amounts to setting the x,y coordinates in all qubits in all gates to 0.
        """
        gates = deepcopy(self.gates)

        name = copy(self.name)

        for i, gate in enumerate(gates):
            gate.pre_wrapped = self.gates[i]
            for qubit in gate.qubits:
                qubit.x = 0
                qubit.y = 0

        circ = Circ(gates)
        circ.name = name

        return circ

    def get_critical_path_length(self):
        """
        Return the number of gates in the longest critical path in the DAG reprensentation of the wrapped circuit.
        """
        wc = self.to_wrapped_circ()
        dag = wc.to_dag()
        lp = nx.dag_longest_path(dag)
        assert (
            type(lp[0][0]) == type(lp[-1][0]) == Qubit
        ), "Init and final nodes of a DAG have to be qubits." + str(type(lp[0][0]))
        return len(lp) - 2  #  We need the number of gates in the longest path.

    def get_gate_dependencies(self):
        """
        Return list l with entries of the form (a,b), indicating that gate a comes before gate b. Only direct dependencies are included.
        """
        wc = self.to_wrapped_circ()
        dag = wc.to_dag()
        deps = []
        for gate in dag:
            if (
                type(gate) == Gate
            ):  # The init and final nodes are qubits, so select only the gates here.
                for pred in dag.predecessors(gate):
                    if type(pred) == Gate:
                        assert (
                            pred.pre_wrapped in self.gates
                        ), "{} not in Basis Circuit".format(pred.pre_wrapped)
                        assert (
                            gate.pre_wrapped in self.gates
                        ), "{} not in Basis Circuit".format(gate.pre_wrapped)

                        dep = (pred.pre_wrapped, gate.pre_wrapped)
                        deps.append(dep)
        return deps

    def rescheduled(self):
        """
        Return a basis circuit with all gates scheduled as early as possible while retaining tilability.
        """
        wbc = self.to_wrapped_circ()
        wbcd = wbc.to_dag()
        wbc = BasisCirc.from_dag(wbcd)
        gates = []
        for gate in wbc.gates:
            gate.qubits = gate.pre_wrapped.qubits
            gates.append(gate)
        bc = BasisCirc(gates, name=self.name)
        return bc

    def merged_gates(self):
        """
        Return copy where  any single- or two-qubit subsequent gates acting on the same qubits are merged. The later gate is merged into the earlier gate. Reschedule after merging.

        If the earlier gate acts on qubits (a,b) and the later gate as well, the new gate gets the name `<later gate name>_<earlier gate name>`. If the earlier gate acts on qubits (a,b) and the later gate on qunits (b,a) than the new gate gets the name `<later gate name>_interchanged_<earlier gate name>`.

        TODO : The function only looks at the current time and the following. Implement using the DAG representation.
        """
        bc = self.rescheduled()
        layers = self.get_layers()
        p = True
        for t in range(len(layers) - 1):
            for gate in layers[t]:
                assert gate.t == t
                for gatep in layers[t + 1]:
                    eq = gatep.qubits == gate.qubits
                    inv_eq = gatep.qubits[::-1] == gate.qubits
                    is_rswap = gatep.name == "rSWAP"
                    if eq or inv_eq:
                        if p == True:
                            print(
                                "NOTE: merging gates and gate names. Also updating the QASM unitary is not yet implemented"
                            )
                            p = False
                        if is_rswap and gate.name == "rSWAP":
                            layers[t + 1].remove(gatep)
                        elif is_rswap:
                            assert gatep.swap == False
                            gate.swap = not gate.swap
                        elif gatep.swap == True:
                            gate.name = f"{gatep.name}_{gate.name}"
                            gate.swap = not gate.swap
                        elif eq:
                            gate.name = f"{gatep.name}_{gate.name}"
                        elif inv_eq:
                            gate.name == f"{gatep.name}_interchanged_{gate.name}"
                        layers[t + 1].remove(gatep)
        gates = [gate for layer in layers for gate in layer]
        bc = BasisCirc(
            gates,
            self.name,
            self.collision_checking,
            self.unique_seeds_checking,
            self.congruent_seeds_checking,
            self.wrapped_circ_checking,
            self.qubits,
        )
        bc = bc.rescheduled()
        return bc

    def to_suzuki_basis_circ(self, r, order=1, merge_gates=True):
        """
        Create basis circuit according to 1st order (order=1) or second order (order=2) suzuki formula with trotter number r. (r repetitions of the cirucit)
        """
        if order == 1:
            assert self.phys_qubits[0] == self.phys_qubits[-1]
            bc = self.repeated(r)
            if merge_gates == True:
                bc = bc.merged_gates()

        elif order == 2:
            print("Creating second order circuit")
            layers = self.get_layers()
            rev_layers = deepcopy(layers)[::-1]
            tmax = self.get_tmax()
            for t, layer in enumerate(rev_layers):
                for gate in layer:
                    gate.t = t + (tmax + 1)
            layers = layers + rev_layers
            for t, layer in enumerate(layers):
                for gate in layer:
                    assert gate.t == t
            gates = [gate for layer in layers for gate in layer]
            bc = BasisCirc(
                gates,
                name=self.name,
                collision_checking=self.collision_checking,
                unique_seeds_checking=self.unique_seeds_checking,
                congruent_seeds_checking=self.congruent_seeds_checking,
                wrapped_circ_checking=self.wrapped_circ_checking,
                qubits=self.qubits,
            )
            bc = bc.repeated(r)
            if merge_gates == True:
                bc = bc.merged_gates()
        else:
            raise ValueError

        if self.phys_qubits != None:
            # Reconstruct phys_qubits
            layers = bc.get_layers()
            init_map = self.phys_qubits[0]
            phys_qubits = deepcopy([init_map])
            for layer in layers:
                map_layer = deepcopy(phys_qubits[-1])
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

                phys_qubits.append(map_layer)

            assert phys_qubits[0] == self.phys_qubits[0]
            assert phys_qubits[-1] == self.phys_qubits[0]
            if isinstance(self, RoutedBasisCirc):
                bc = RoutedBasisCirc(
                    bc.gates,
                    phys_qubits,
                    name=f"self.name_Suzuki_r={r}_order={order}",
                    collision_checking=self.collision_checking,
                    unique_seeds_checking=self.unique_seeds_checking,
                    congruent_seeds_checking=self.congruent_seeds_checking,
                    wrapped_circ_checking=self.wrapped_circ_checking,
                    qubits=self.qubits,
                )
            else:
                bc = BasisCirc(
                    bc.gates,
                    name=f"self.name_Suzuki_r={r}_order={order}",
                    collision_checking=self.collision_checking,
                    unique_seeds_checking=self.unique_seeds_checking,
                    congruent_seeds_checking=self.congruent_seeds_checking,
                    wrapped_circ_checking=self.wrapped_circ_checking,
                    qubits=self.qubits,
                )
                bc.phys_qubits = phys_qubits

        return bc


class RoutedBasisCirc(BasisCirc):
    """
    A basis circ with no congruent seeds checking by default. Also, get_patch will automaticaly return the completed patch, and get_gates will return a completed patch of 1 by 1 basis circuits. Also, providing a qubit map phys_qubits is mandatory.

    TODO : Make sure all calls to BasisCirc are to RoutedBasisCirc in the apropriate places.
    """

    def __init__(
        self,
        gates,
        phys_qubits,
        name=None,
        collision_checking=True,
        unique_seeds_checking=True,
        congruent_seeds_checking=False,  # Different from BasisCirc
        wrapped_circ_checking=True,
        qubits=None,
    ):
        super().__init__(
            gates,
            name=name,
            collision_checking=collision_checking,
            unique_seeds_checking=unique_seeds_checking,
            congruent_seeds_checking=congruent_seeds_checking,
            wrapped_circ_checking=wrapped_circ_checking,
            qubits=qubits,
        )
        self.phys_qubits = phys_qubits

    def get_gates(self, n, m, nc=1):
        p = self.get_completed_patch(n, m, phys_qubits=self.phys_qubits, nc=nc)
        return p.gates

    def get_patch(self, n, m, nc=1):
        return self.get_completed_patch(n, m, phys_qubits=self.phys_qubits, nc=nc)


class Edge(Gate):
    """
    Just as a gate but used to describe device connectivity.
    """

    def __init__(self, qubits, name="EDGE", t=None):
        super().__init__(name, qubits, t=t)
        assert self.n == 2

    def __eq__(self, other):
        return (
            type(other) == type(self)
            and self.name == other.name
            and (self.qubits == other.qubits or self.qubits == other.qubits[::-1])
        )  # Ignore order of qubits

    def to_tuple(self):
        return (self.qubits[0].to_tuple(), self.qubits[1].to_tuple())

    def __hash__(self):
        return hash(
            (self.name, tuple(sorted(self.to_tuple())))
        )  # Hashing ignores qubit order. (I.e. edges are undirected)


class RSwap(Edge):
    """
    Routing SWAP. Just as Edge, but now has an extra attribute `swap` that says whether the SWAP is 'on' (whether it is actually performed or not). This attrubite can be a z3.Bool. Hashing ignores qubit order self.t and self.on. Same for "==".
    """

    def __init__(self, qubits, t, on):
        super().__init__(qubits, name="rSWAP", t=t)
        self.on = on


class BasisGraph(BasisCirc):
    """
    A BasisCirc, but used to describe device connectivity. It has Edges instead of Gates.

    `edges` is a set of `Edge`s. All seed numbers s together must form a congruent set of integers. Every seed in the unit cell must have an edge attached to it.
    """

    def __init__(self, edges, name=None):
        self.edges = set(edges)
        super().__init__(
            list(self.edges),
            name=name,
            collision_checking=False,
            unique_seeds_checking=False,
            congruent_seeds_checking=True,
            wrapped_circ_checking=False,
        )
        assert (
            self.no_lone_seeds()
        ), "Please redefine the basis graph so that every seed (in the base unit cell) is attached to an edge (this is WLOG)."

    def no_lone_seeds(self):
        """
        Return `True` if all seed numbers appear in the unit cell.
        """
        smax = self.get_smax()
        ucqs = self.get_seeds()
        seed_nums_in_uc = {q.s for q in ucqs}
        return sorted(list(seed_nums_in_uc)) == list(range(smax + 1))

    def to_basis_circ(self, random=False, long_gate_names=True):
        """
        Return BasisCirc where one two-qubit gate with name exp_H_ij is placed along every edge ij of self. Every gate occurs at a new time step. Gates are placed in the order of `list(basis_graph.edges)`. If (v,v') is an edge of the basis raph, the control of the two-qubit gate is v, and v' is the target.
        If `random==True`, each gate name will be gate_name(rand) with rand a number between -pi and pi chosen at random for every gate.
        """
        gates = []
        for t, edge in enumerate(self.edges):
            q0, q1 = deepcopy(edge.qubits)
            if q0.x < 0:
                q0.x = "min" + str(abs(q0.x))
            if q0.y < 0:
                q0.y = "min" + str(abs(q0.y))
            assert q0.s >= 0
            if q1.x < 0:
                q1.x = "min" + str(abs(q1.x))
            if q1.y < 0:
                q1.y = "min" + str(abs(q1.y))
            assert q1.s >= 0

            if long_gate_names == True:
                gn = f"exp_H_{q0.x}_{q0.y}_{q0.s}__{q1.x}_{q1.y}_{q1.s}"
                if random == True:
                    rn = np.random.rand() * 2 * np.pi - np.pi
                    gn += f"({rn})"
            else:
                gn = "G"
            gate = Gate(gn, edge.qubits, t)
            gates.append(gate)

        bc = BasisCirc(gates, name=self.name)
        return bc

    @classmethod
    def from_json(cls, name, path=files("quantile").joinpath("basis_graphs.json")):
        with open(path, "r") as f:
            db = json.load(f)

        assert name in db, "{} not in {}".format(name, path)
        edges = db[name]["edges"]
        edges = map(ast.literal_eval, edges)
        edges = [Edge((Qubit(*edge[0]), Qubit(*edge[1]))) for edge in edges]

        return cls(edges, name=name)

    def get_plot_graph(self, node_labels=False):
        """
        Return a patch of 3x3 basis graphs as a networkx graph, with the middle basis graph bold.
        """

        def cartesian_to_string(c):
            spos = f"{c[0]}.0,{c[1]}.0"
            return spos

        penwidth = 5
        wide_penwidth = 16
        width = 0.1
        headclip = "false"
        tailclip = "false"
        base_color = "#80808099"  # Translucent Gray
        fontsize = 10

        p = self.get_open_patch(5, 5)
        tp = p.translated(-2, -2)

        edges = [edge.to_tuple() for edge in tp.gates]
        g = nx.Graph(edges)

        nx.set_node_attributes(g, width, "width")
        nx.set_edge_attributes(g, headclip, "headclip")
        nx.set_edge_attributes(g, tailclip, "tailclip")
        nx.set_edge_attributes(g, penwidth, "penwidth")
        nx.set_edge_attributes(g, base_color, "color")
        nx.set_node_attributes(g, fontsize, "fontsize")
        nx.set_node_attributes(g, "plaintext", "shape")
        if node_labels == True:
            mapping = {node: {"label": node[:3]} for node in g}
            nx.set_node_attributes(g, mapping)
        else:
            nx.set_node_attributes(g, "point", "shape")

        for edge in tp.gates:
            if edge in self.edges:
                edge = edge.to_tuple()
                g[edge[0]][edge[1]]["penwidth"] = wide_penwidth

        assert all(
            type(node[0]) == type(node[1]) == int for node in g
        ), "First two entries of a node must be ints, representing the x,y coordinate of the cell the node is in."
        mapping = {node: {"pos": cartesian_to_string(node[:2])} for node in g}
        nx.set_node_attributes(g, mapping)

        return g

    def save_plot(self, folder="../basis_graph_plots", node_labels=False):
        g = self.get_plot_graph(node_labels=node_labels)
        g = nx.nx_agraph.to_agraph(g)
        import warnings

        warnings.filterwarnings("ignore", category=RuntimeWarning, module="pygraphviz")
        g.draw(
            folder + "/" + self.name + ".pdf",
            prog="neato",
            args="-Goutputorder=edgesfirst",
        )


def load_solution_database(fname="solutions.pkl"):
    import pickle

    if os.path.isfile(fname):
        with open(fname, "rb") as f:
            db = pickle.load(f)
    else:
        db = {}

    return db


def plot_dag(c, name):
    """
    Plot the DAG circuit c (networkx.MultiDiGraph), and save under `name`. File extension is set by extension of name.
    """
    g = deepcopy(c)
    for _, _, data in g.edges(data=True):
        data["label"] = data["qubit"]

    g = nx.nx_agraph.to_agraph(g)

    g.draw(
        name,
        prog="dot",
        args="-Grankdir=BT",
    )
