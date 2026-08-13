"""
End-to-end pipeline for the obstructing n-friend pairs in
Data/results/obstructed_pairs.csv.

For every row (or a --pairs subset) this script:
  1. Reads the census-knot PD code and the n-friend PD code.
  2. Calls forms_special_NRBG_link(n, blue_ex, green_ex) from code/n_rbg.py to
     find an n-special RBG link over the pair.
  3. On success, writes into Data/results/example_knots/<pair>/:
       - K_census.pd.txt          PD code of the census knot K
       - K_friend.pd.txt          PD code of the n-friend K'
       - RBG_link.pd.txt          PD code of the 3-component RBG link
       - README.md                id_num, n, invariants, source rows, framings
       - K_census.png / .pdf      orthogonal projection (single-component: black)
       - K_friend.png / .pdf
       - RBG_link.png / .pdf      (R, B, G colors match the component order)

Folders whose RBG_link.pd.txt already exists are skipped unless --force is
passed, so re-running resumes cleanly after an interruption.

Requires Sage (for n_rbg.py's `vector`/`matrix`), snappy, plink, and
`pdftoppm` (from poppler) on PATH.

Run from the project root:

    conda activate sage
    python code/build_example_knots.py                    # all obstructing pairs
    python code/build_example_knots.py --pairs 20,166,168 # subset by id_num
    python code/build_example_knots.py --force            # redo even if present
    python code/build_example_knots.py --skip-render      # PD/README only
"""

import argparse
import ast
import csv
import os
import sys
import traceback

import snappy

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

# n_rbg.py depends on Sage (it uses unqualified `vector`, `matrix`, `Slope`),
# and it is not written as an importable module — it declares free-floating
# helpers relying on Sage globals. So we exec it in this module's globals
# rather than `import`-ing it, and we pull the Sage builtins in first.
from sage.all import matrix, vector  # noqa: F401 (used inside n_rbg.py)
try:
    from sage.all import Slope  # noqa: F401
except ImportError:
    Slope = object  # only used as a type hint in n_rbg.py

with open(os.path.join(_HERE, "n_rbg.py")) as _f:
    exec(compile(_f.read(), "n_rbg.py", "exec"), globals())

# Local renderer (also import lazily so --skip-render still works if plink
# is unavailable).
try:
    from render_example_knots import render_pd_file  # noqa: E402
    _RENDER_IMPORT_ERROR = None
except Exception as _e:  # pragma: no cover
    render_pd_file = None
    _RENDER_IMPORT_ERROR = _e

MASTER_CSV = os.path.join(_ROOT, "Data", "n_friends(master version).csv")
OBSTRUCTED_CSV = os.path.join(_ROOT, "Data", "results", "obstructed_pairs.csv")
OUT_DIR = os.path.join(_ROOT, "Data", "results", "example_knots")


def _parse_pd(s):
    return ast.literal_eval(s)


def _load_all_rows():
    with open(OBSTRUCTED_CSV) as f:
        return list(csv.DictReader(f))


def _folder_name(row):
    return f"id{row['id_num']}_n{row['n']}_{row['original_knot_name']}"


def _write_pd(link, path):
    with open(path, "w") as f:
        f.write(repr(link.PD_code()))
        f.write("\n")


def _write_readme(path, row, framings):
    id_num = row["id_num"]
    n = row["n"]
    census = row["original_knot_name"]
    lines = [
        f"# Example pair: id_num {id_num} (n = {n})",
        "",
        f"- Source: `Data/n_friends(master version).csv` row `id_num = {id_num}`",
        f"- Companion: `Data/results/obstructed_pairs.csv`",
        f"- Census knot K: `{census}` "
        f"(index {row['n_friend_index']} in `plausibly_unknown.csv`, "
        f"{row['original_knot_crossings']} crossings, "
        f"τ = {row['original_knot_tau']}, s₀ = {row['original_knot_s_0']})",
        f"- n-friend K': {row['num_crossings']} crossings; "
        f"τ = {row['tau']}, ν = {row['nu']}, "
        f"s = {row['s']}, s₂ = {row['s_2']}, s₃ = {row['s_3']}",
        f"- Obstruction: {row['obstructing_invariants']}",
        f"- RBG framings (r, b, g): {framings}",
        "",
        "## Files",
        "- `K_census.{pd.txt,png,pdf}` — census knot K (single component, black).",
        "- `K_friend.{pd.txt,png,pdf}` — its n-friend K' (single component, black).",
        "- `RBG_link.{pd.txt,png,pdf}` — n-special RBG link (R = red, B = blue, G = green).",
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _folder_is_complete(folder, want_png):
    needed = ["K_census.pd.txt", "K_friend.pd.txt", "RBG_link.pd.txt", "README.md"]
    if want_png:
        needed += ["K_census.png", "K_friend.png", "RBG_link.png"]
    return all(os.path.exists(os.path.join(folder, name)) for name in needed)


def _render_folder(folder):
    if render_pd_file is None:
        print(f"    (skip render — import error: {_RENDER_IMPORT_ERROR})")
        return False
    ok = True
    for name in ("K_census.pd.txt", "K_friend.pd.txt", "RBG_link.pd.txt"):
        try:
            render_pd_file(os.path.join(folder, name))
        except Exception as e:  # pragma: no cover
            print(f"    ✗ render {name}: {type(e).__name__}: {e}")
            ok = False
    return ok


def _process_row(row, force, do_render):
    folder = os.path.join(OUT_DIR, _folder_name(row))
    rbg_pd_path = os.path.join(folder, "RBG_link.pd.txt")

    header = f"=== id_num {row['id_num']} (n = {row['n']}) over {row['original_knot_name']} ==="
    print(f"\n{header}")

    already_searched = os.path.exists(rbg_pd_path)
    if already_searched and not force:
        if _folder_is_complete(folder, want_png=do_render):
            print("  ↷ already complete, skipping")
            return "skipped"
        if do_render:
            print("  ↷ PD codes present, (re)rendering PNGs only")
            _render_folder(folder)
            return "rendered"

    n = int(row["n"])
    census_pd = _parse_pd(row["n_friend_PD_code"])
    friend_pd = _parse_pd(row["knot_PD_code"])
    K_census = snappy.Link(census_pd)
    K_friend = snappy.Link(friend_pd)

    print(f"  K crossings  = {len(K_census.crossings)}")
    print(f"  K' crossings = {len(K_friend.crossings)}")
    print("  Searching for n-special RBG link …")

    rbg = forms_special_NRBG_link(n, K_census.exterior(), K_friend.exterior())  # noqa: F821
    if rbg is None:
        print("  ✗ no n-special RBG link found")
        return "failed"

    os.makedirs(folder, exist_ok=True)
    _write_pd(K_census, os.path.join(folder, "K_census.pd.txt"))
    _write_pd(K_friend, os.path.join(folder, "K_friend.pd.txt"))
    _write_pd(rbg.link, os.path.join(folder, "RBG_link.pd.txt"))
    _write_readme(os.path.join(folder, "README.md"), row, rbg.framings)
    print(f"  ✓ wrote {os.path.relpath(folder, _ROOT)}")

    if do_render:
        _render_folder(folder)
    return "found"


def _parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--pairs",
        default="all",
        help="Comma-separated id_num list, or 'all' (default). "
        "Example: --pairs 20,166,168",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Re-run the RBG search even if an output folder already exists.",
    )
    p.add_argument(
        "--skip-render",
        action="store_true",
        help="Skip PDF/PNG rendering; write only PD codes and README.",
    )
    return p.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv or sys.argv[1:])
    os.makedirs(OUT_DIR, exist_ok=True)

    all_rows = _load_all_rows()
    if args.pairs.strip().lower() == "all":
        rows = all_rows
    else:
        wanted = {p.strip() for p in args.pairs.split(",") if p.strip()}
        rows = [r for r in all_rows if r["id_num"] in wanted]
        missing = wanted - {r["id_num"] for r in rows}
        if missing:
            print(f"warning: id_num(s) not found in obstructed_pairs.csv: "
                  f"{sorted(missing)}")

    print(f"Processing {len(rows)} pair(s). Output → {OUT_DIR}")

    results = {}
    for row in rows:
        try:
            results[row["id_num"]] = _process_row(
                row, force=args.force, do_render=not args.skip_render
            )
        except KeyboardInterrupt:
            print("\n^C — aborting; partial output preserved.")
            results[row["id_num"]] = "interrupted"
            break
        except Exception as e:
            print(f"  ✗ id_num {row['id_num']} raised "
                  f"{type(e).__name__}: {e}")
            traceback.print_exc()
            results[row["id_num"]] = "error"

    print("\n=== Summary ===")
    by_status = {}
    for id_num, status in results.items():
        by_status.setdefault(status, []).append(id_num)
    for status in ("found", "skipped", "rendered", "failed", "error", "interrupted"):
        ids = by_status.get(status, [])
        if ids:
            print(f"  {status:<11} {len(ids):>3}  {', '.join(ids)}")
    return results


if __name__ == "__main__":
    main()
