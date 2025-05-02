#!/us4,4r/bin/env python3

import sys
import re

sys.path.insert(1, "../")
import quantile as qt
from tabulate import tabulate
from itertools import groupby
import json


db = qt.load_solution_database("solutions.pkl")


def partition_dict_list(dict_list, exclude_keys):
    # Helper function to create a "comparison key" by excluding the specified key
    def without_keys(d, exclude_keys):
        d = {k: v for k, v in d.items() if k not in exclude_keys}
        return list(d.values())

    # Sort the list to prepare for grouping
    dict_list_sorted = sorted(dict_list, key=lambda d: without_keys(d, exclude_keys))

    # Group by the dictionaries excluding the exclude keys.  Sort the groups by sort_key
    partitions = []
    for _, group in groupby(
        dict_list_sorted, key=lambda d: without_keys(d, exclude_keys)
    ):
        partitions.append(list(group))  # Convert group iterator to list

    return partitions


def split(s):
    index = s.find("(")

    if index == -1:  # If '(' is not found, return the original string and an empty part
        return s, ""

    left = s[:index]
    right = s[index:]

    return left, right


def calc_lines(group_name, transpiler_name=None):
    lines = []
    if group_name == "fermi-hubbard":
        group_name = "Fermi-Hubbard"
    for t in db[group_name]:
        if (transpiler_name == None or t.name == transpiler_name) and t.solution[
            "solved"
        ] == True:
            print(t.basis_circuit.get_active_qubits())
            bc = t.basis_circ
            bc = bc.resized(3, 3)[0]  # Assumes a reach of max 1.
            if group_name == "two-qudit":
                wbc = bc.to_wrapped_circ()
                min_depth = (
                    wbc.get_max_qubit_activity()
                )  # For two-qubit, assumes graph is Vizing class I
            elif group_name == "kogut-susskind":
                min_depth = 210  # Always same depth anyway
            else:
                min_depth = bc.get_critical_path_length()
            depth = t.solution["depth"]
            depth_overhead = depth - min_depth
            depth_overhead_pct = depth_overhead / min_depth * 100
            bg_num_seeds = t.basis_graph.get_smax() + 1
            qubit_overhead = bg_num_seeds - (t.basis_circ.get_smax() + 1)
            print(t.basis_circ.name)
            bc_name, bc_size = split(t.basis_circ.name)
            bg_name, bg_size = split(t.basis_graph.name)
            rbc = t.solution["routed_basis_circ"]
            nkd_swaps = len(
                list(g for g in rbc.gates if g.name == "rSWAP")
            )  # Number of naked swaps per routed basis circ.

            try:
                fixed_naked_swaps = t.fixed_naked_swaps  # For backward comp.
            except:
                fixed_naked_swaps = False

            lines.append(
                {
                    "basis circ": bc_name,
                    "bc size": bc_size,
                    "basis graph": bg_name,
                    "bg size": bg_size,
                    "min_depth": min_depth,
                    "depth": depth,
                    "d. overh.": depth_overhead,
                    "d. overh. (%)": depth_overhead_pct,
                    "merge swaps": t.merge_swaps,
                    "cyclic": t.cyclic,
                    "q. overh.": qubit_overhead,
                    "nkd swaps": nkd_swaps,
                    "bg num seeds": bg_num_seeds,
                    "minimize swaps": t.minimize_swaps,
                    "fixed_naked_swaps": fixed_naked_swaps,
                    "wall clock": t.solution["wall_clock"],
                    "slice depth": t.slice_depth,
                }
            )
    for line in lines:
        print(line)
    return lines


def filtered_lines(lines):
    def squareness(line):  # lower is squarer
        if len(line["bc size"]) == 5 and len(line["bg size"]) == 5:
            n, m = [int(line["bc size"][1]), int(line["bc size"][3])]
            nn, mm = [int(line["bg size"][1]), int(line["bg size"][3])]
            sqns = abs(n - m) + abs(nn - mm)
        else:
            sqns = 0
        return sqns

    # Do not consider heavy-hex and square octagon, since they are monomorphisms of square lattice.
    exclude = ["heavy-hex", "square-octagon"]
    lines = [line for line in lines if line["basis circ"] not in exclude]

    # Keep only the size with the least depth overhead,
    # and of those, with the least number of naked swaps,
    # and of those, with the least number of seeds in the basis graph,
    # and of those, the one where minimize swaps was on,
    # and of those, the one with the highest squareness,
    # and of those, the one that was fastest.
    ignore = [  # What is allowed to be different per partition from which the best is selected.
        "bc size",
        "bg size",
        "depth",
        "min_depth",
        "d. overh.",
        "d. overh. (%)",
        "q. overh.",
        "bg num seeds",
        "minimize swaps",
        "nkd swaps",
        "fixed_naked_swaps",
        "wall clock",
    ]

    _lines = []
    partitions = partition_dict_list(lines, ignore)
    for part in partitions:
        _part = sorted(
            part,
            key=lambda x: [
                x["d. overh."],
                x["nkd swaps"],
                x["bg num seeds"],
                not x["minimize swaps"],
                squareness(x),
                x["wall clock"],
            ],
        )
        _lines.append(_part[0])

    lines = _lines

    # Sort the lines
    def key(line):
        return [
            line["d. overh."],
            line["nkd swaps"],
            line["q. overh."],
            line["merge swaps"],
            not line["cyclic"],
        ]

    lines.sort(key=key)
    return lines


def to_dhms(time):
    def pad(des):
        if len(des) == 1:
            des = "0" + des
        return des

    day = time // (24 * 3600)
    time = time % (24 * 3600)
    hour = time // 3600
    time %= 3600
    minutes = time // 60
    time %= 60
    seconds = time
    f = [day, hour, minutes, seconds]
    f = [pad(str(des)) for des in f]
    string = f"{f[0]}:{f[1]}:{f[2]}:{f[3]}"
    string = string[1:]
    return string


def latex_table(group_name, json_file=None):
    with open(json_file, "r") as f:
        lines = json.load(f)

    lines = filtered_lines(lines)

    for line in lines:
        do = line["d. overh."]
        dopct = round(line["d. overh. (%)"])
        line[
            "do (do %)"
        ] = f"{do} ({dopct} %)"  # Include relative overhead in overhead field.
        line["wall clock"] = round(line["wall clock"])
        line["wall clock"] = to_dhms(line["wall clock"])

        if (
            line["minimize swaps"] == False
        ):  # Mark the lines where swaps were not minimized.
            line["nkd swaps"] = f"{line['nkd swaps']}*"
        if (
            type(line["fixed_naked_swaps"]) == int and line["merge swaps"] == True
        ):  # Mark the lines run at a fixed number of naked swaps
            assert line["nkd swaps"] == line["fixed_naked_swaps"]
            line["wall clock"] = f"{line['wall clock']}*"
        line["wall clock"] = "texttt{" + line["wall clock"] + "}"

    keys_to_include = [
        "basis circ",
        "bc size",
        "basis graph",
        "bg size",
        "do (do %)",
        "nkd swaps",
        "q. overh.",
        "merge swaps",
        "cyclic",
        "wall clock",
    ]
    custom_headers = [
        "Basis circuit",
        "Size",
        "Basis graph",
        "Size",
        "depth\noverhead",
        "naked\nswaps",
        "qubit\noverhead",
        "merge\nswaps",
        "cyclic",
        "wall clock",
    ]

    if group_name == "kogut-susskind":
        keys_to_include.insert(-1, "slice depth")
        custom_headers.insert(-1, "slice depth")

    # Extract the relevant columns for the table
    filtered_data = [[l[key] for key in keys_to_include] for l in lines]

    table = tabulate(filtered_data, headers=custom_headers, tablefmt="latex_longtable")
    table = table.replace(r"\{", r"{")
    table = table.replace(r"\}", r"}")
    table = table.replace(r"texttt", r"\texttt")
    return table


def save_data(group_name):
    print(group_name)
    lines = calc_lines(group_name)
    with open(f"{group_name}/data.json", "w") as f:
        json.dump(lines, f, indent=4)


def make_table(group_name):
    table = latex_table(group_name, f"{group_name}/data.json")
    pattern = r"\\endhead\s*(.*?)\s*\\hline\n\\end\{longtable\}"
    match = re.search(pattern, table, re.DOTALL)
    assert match
    table_body = match.group(1).strip()  # Remove any leading/trailing whitespace

    with open(f"{group_name}/{group_name}_table_body.tex", "w") as f:
        f.write(table_body)


# name = "kogut-susskind"
# name = "rule54"
# name='two-qubit'
# name = "rokhsar-kivelson"
name = "fermi-hubbard"
save_data(name)
make_table(name)
