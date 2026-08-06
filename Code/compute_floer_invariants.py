"""Computes Tau and Nu invariants from planar diagram codes."""

from __future__ import annotations

import argparse
import ast
import csv
import io
import os
import re
import shutil
import subprocess
import tempfile
from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd
import math
import time

import snappy

from invariant_io import resolve_data_path, update_row, value_is_missing

try:
    import caffeine
except ImportError:
    caffeine = None

# Override this with N_FRIENDS_DATA or the data_path function argument.
out_file_path = str(resolve_data_path())

data_types = {
    "id_num": int,
    "num_crossings": int,
    "volume": float,
    "n": int,
    "verification": bool,
    "s": str,
    "s_2": str,
    "s_3": str,
    "nu": str,
    "tau": str,
    "obstructs": bool,
    "n_friend_name": str,
    "n_friend_index": int,
    "knot_PD_code": str,
    "n_friend_PD_code": str
    }

def compute_floer_homology(
    id_num: int,
    max_crossings: int = 100,
    recompute: bool = False,
    data_path=None,
):
    path = resolve_data_path(data_path)
    data = pd.read_csv(path, dtype=data_types)
    matches = data.loc[data["id_num"] == int(id_num)]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one row with id_num={id_num}, found {len(matches)}.")
    row = matches.iloc[0]

    #print(row)
    num_crossings = int(row["num_crossings"])

    if num_crossings > max_crossings:
        #print("Too many crossings")
        return
    
    pd_code = ast.literal_eval(row["knot_PD_code"])
    n = int(row["n"])
    if not (value_is_missing(row["nu"]) or recompute):
        return

    K=snappy.Link(pd_code)

    #Computes Knot Floer Homology
    print(f"Finding Knot Floer Homology for knot {id_num}")

    start_time = time.perf_counter()


    homology = K.knot_floer_homology()

        
    print("nu=" + str(homology["nu"]))
    print("tau=" + str(homology["tau"]))

    elapsed = time.perf_counter() - start_time

    hours, remainder = divmod(elapsed, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"Runtime: {int(hours)}h {int(minutes)}m {seconds:.2f}s")

    nu = homology["nu"]
    tau = homology["tau"]
    
    updates = {"nu": str(nu), "tau": str(tau)}

    #This flips the sign of the inequalities we check when n is negative
    #Note: this doesn't change the sign of what's written to the data file
    if (n < 0):
        nu = -nu
        tau = -tau
        n = -n

    obstructs = False
    #Checks Tau obstruction
    if (2*tau > n-math.sqrt(n)):
        obstructs = True

    #Checks Nu obstruction (double check that this is correct)
    if (2*nu > n-math.sqrt(n)):
        obstructs = True
    
    if obstructs:
        updates["obstructs"] = True

    update_row(id_num, updates, data_path=path)
    return updates

def compute_floer_invariants_specified(indices: list[int], max_crossings: int = 100, data_path=None):
    if caffeine:
        caffeine.on()
    for i in indices:
        compute_floer_homology(i, max_crossings, data_path=data_path)
    if caffeine:
        caffeine.off()

def compute_floer_invariants_table(start: int, end: int, max_crossings: int = 100, recompute: bool =False, data_path=None):
    if caffeine:
        caffeine.on()
    for i in range(start,end + 1):
        compute_floer_homology(i, max_crossings, recompute, data_path=data_path)
    if caffeine:
        caffeine.off()
