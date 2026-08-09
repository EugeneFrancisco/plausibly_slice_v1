"""
Build example RBG-link/n-friend folders under
preliminary_results/example_knots/.

Run in Sage from the project root:

    sage: load('Code/build_example_knots.py')
    sage: main()

For each target row in preliminary_results/obstructed_pairs.csv this script:
  1. Reads the census-knot PD code and the n-friend PD code.
  2. Calls forms_special_NRBG_link(n, blue_ex, green_ex) from Code/n_rbg.py to
     find an n-special RBG link over the pair.
  3. On success, writes into a per-pair folder:
       - K_census.pd.txt          PD code of the census knot K
       - K_friend.pd.txt          PD code of the n-friend K'
       - RBG_link.pd.txt          PD code of the 3-component RBG link
       - K_census.lnk / .lnk      best-effort snappy Link projections
       - K_friend.lnk / .lnk
       - RBG_link.lnk
       - README.md                id_num, n, invariants, source rows, framings

The .lnk save is best-effort: it needs plink's Tk editor. If it fails
(e.g. running headless), only the PD-code text files are written and the
README notes so.
"""

import ast
import csv
import os
import sys
import traceback

import snappy

# n_rbg.py depends on Sage (it uses `vector`, `matrix`, `Slope`, ...), and it
# is not written as an importable module — it declares free-floating helpers
# and Sage builtins. So we exec it in this module's globals rather than
# `import`-ing it, and we pull the Sage builtins in first so its unqualified
# uses of `matrix` / `vector` / `Slope` resolve.
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from sage.all import matrix, vector  # noqa: F401 (used inside n_rbg.py)
try:
    from sage.all import Slope  # noqa: F401
except ImportError:
    Slope = object  # only used as a type hint in n_rbg.py

with open(os.path.join(_HERE, "n_rbg.py")) as _f:
    exec(compile(_f.read(), "n_rbg.py", "exec"), globals())

MASTER_CSV = os.path.join(_ROOT, "data", "n_friends(master version).csv")
OBSTRUCTED_CSV = os.path.join(_ROOT, "preliminary_results", "obstructed_pairs.csv")
OUT_DIR = os.path.join(_ROOT, "preliminary_results", "example_knots")

# id_nums picked from preliminary_results/obstructed_n_friends.md as the three
# smallest / most-tractable obstructing pairs. All have n > 0, so no mirroring
# of the PD codes is required before feeding them to snappy.
TARGET_IDS = ["20", "166", "168"]


def _parse_pd(s):
    return ast.literal_eval(s)


def _load_row(id_num):
    with open(OBSTRUCTED_CSV) as f:
        for row in csv.DictReader(f):
            if row["id_num"] == id_num:
                return row
    raise KeyError(f"id_num {id_num} not in obstructed_pairs.csv")


def _folder_name(row):
    return f"id{row['id_num']}_n{row['n']}_{row['original_knot_name']}"


def _write_pd(link, path):
    with open(path, "w") as f:
        f.write(repr(link.PD_code()))
        f.write("\n")


def _try_save_lnk(link, path):
    """Attempt to save a plink .lnk projection. Return True on success."""
    try:
        editor = link.view()
        editor.save_file(path)
        try:
            editor.done()
        except Exception:
            pass
        return True
    except Exception:
        try:
            link.save_file(path)
            return True
        except Exception:
            return False


def _write_readme(path, row, framings, saved_lnk):
    n = row["n"]
    id_num = row["id_num"]
    census = row["original_knot_name"]
    obstruction = row["obstructing_invariants"]
    lines = [
        f"# Example pair: id_num {id_num} (n = {n})",
        "",
        f"- Source CSV: `data/n_friends(master version).csv` row `id_num = {id_num}`",
        f"- Companion row: `preliminary_results/obstructed_pairs.csv`",
        f"- Census knot K: `{census}` (index {row['n_friend_index']} in `plausibly_unknown.csv`)",
        f"- Original K crossings: {row['original_knot_crossings']}; "
        f"K τ = {row['original_knot_tau']}; K s₀ = {row['original_knot_s_0']}",
        f"- n-friend K' crossings: {row['num_crossings']}",
        f"- K' invariants: τ = {row['tau']}, ν = {row['nu']}, "
        f"s = {row['s']}, s₂ = {row['s_2']}, s₃ = {row['s_3']}",
        f"- Obstruction: {obstruction}",
        f"- RBG framings (r, b, g): {framings}",
        "",
        "## Files",
        "- `K_census.pd.txt` — PD code of K (the census knot).",
        "- `K_friend.pd.txt` — PD code of K' (its n-friend).",
        "- `RBG_link.pd.txt` — PD code of the n-special RBG link found.",
    ]
    if saved_lnk:
        lines.append("- `*.lnk` — SnapPy plink projections for the same links.")
    else:
        lines.append(
            "- `.lnk` projections were not written (headless / plink unavailable)."
            " Re-open the PD codes with `snappy.Link(pd).view().save_file(...)`"
            " in an interactive session to generate them."
        )
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _try_pair(row):
    n = int(row["n"])
    census_pd = _parse_pd(row["n_friend_PD_code"])
    friend_pd = _parse_pd(row["knot_PD_code"])
    K_census = snappy.Link(census_pd)
    K_friend = snappy.Link(friend_pd)

    blue_ex = K_census.exterior()
    green_ex = K_friend.exterior()

    print(f"\n=== id_num {row['id_num']} (n = {n}) over {row['original_knot_name']} ===")
    print(f"  K crossings  = {len(K_census.crossings)}")
    print(f"  K' crossings = {len(K_friend.crossings)}")
    print("  Searching for n-special RBG link …")

    rbg = forms_special_NRBG_link(n, blue_ex, green_ex)  # noqa: F821
    if rbg is None:
        print("  ✗ no n-special RBG link found for this pair.")
        return False

    folder = os.path.join(OUT_DIR, _folder_name(row))
    os.makedirs(folder, exist_ok=True)

    _write_pd(K_census, os.path.join(folder, "K_census.pd.txt"))
    _write_pd(K_friend, os.path.join(folder, "K_friend.pd.txt"))
    _write_pd(rbg.link, os.path.join(folder, "RBG_link.pd.txt"))

    saved_lnk = True
    saved_lnk &= _try_save_lnk(K_census, os.path.join(folder, "K_census.lnk"))
    saved_lnk &= _try_save_lnk(K_friend, os.path.join(folder, "K_friend.lnk"))
    saved_lnk &= _try_save_lnk(rbg.link, os.path.join(folder, "RBG_link.lnk"))

    _write_readme(
        os.path.join(folder, "README.md"),
        row,
        rbg.framings,
        saved_lnk,
    )
    print(f"  ✓ wrote {folder}")
    return True


def main(ids=None):
    os.makedirs(OUT_DIR, exist_ok=True)
    ids = ids or TARGET_IDS
    results = {}
    for id_num in ids:
        row = _load_row(id_num)
        try:
            results[id_num] = _try_pair(row)
        except Exception as e:
            print(f"  ✗ id_num {id_num} raised {type(e).__name__}: {e}")
            traceback.print_exc()
            results[id_num] = False

    print("\n=== Summary ===")
    for id_num, ok in results.items():
        mark = "✓" if ok else "✗"
        print(f"  {mark} id_num {id_num}")
    return results


if __name__ == "__main__":
    main()
