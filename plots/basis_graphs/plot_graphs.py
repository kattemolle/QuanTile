#!/usr/bin/env python3


import sys

sys.path.insert(1, "../")
import quantile as qt
import json


with open("../quantile/basis_graphs.json", "r") as f:
    db = json.load(f)

for name in db:
    print(name)
    bg = qt.BasisGraph.from_json(name, path="../quantile/basis_graphs.json")
    bg.save_plot(node_labels=True)
