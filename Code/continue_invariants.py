"""Resume missing n-friends invariant calculations from the command line."""

from __future__ import annotations

import argparse
import sys
import traceback

import pandas as pd

from compute_floer_invariants import FLOER_COLUMNS, compute_floer_homology
from compute_s_invariant import KNOTJOB_JAR, S_COLUMNS, compute_s_invariants
from invariant_io import resolve_data_path, value_is_missing


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("floer", "s"))
    parser.add_argument("--data", default=None, help="n-friends CSV (defaults to master)")
    parser.add_argument("--start", type=int, default=1, help="first id_num to inspect")
    parser.add_argument("--end", type=int, default=None, help="last id_num to inspect")
    parser.add_argument("--max-crossings", type=int, default=100)
    parser.add_argument("--timeout", type=float, default=3600, help="per-characteristic s timeout")
    parser.add_argument("--jar", default=str(KNOTJOB_JAR), help="KnotJob.jar path")
    parser.add_argument(
        "--order",
        choices=("id", "crossings"),
        default="id",
        help="work through candidates by id_num, or fewest crossings first",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="report what would be computed, then exit without computing",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    path = resolve_data_path(args.data)
    data = pd.read_csv(path)
    end = args.end if args.end is not None else int(data["id_num"].max())
    candidates = data.loc[
        data["id_num"].between(args.start, end)
        & data["num_crossings"].le(args.max_crossings)
    ]
    columns = FLOER_COLUMNS if args.kind == "floer" else S_COLUMNS
    total_in_range = len(candidates)
    candidates = candidates.loc[
        candidates[list(columns)].apply(
            lambda row: any(value_is_missing(value) for value in row), axis=1
        )
    ]

    # Cheapest-first is worth a lot here: the cost of these invariants climbs
    # steeply with crossing number, so ordering by id can park the worker on one
    # enormous knot while hundreds of small ones remain untouched.
    if args.order == "crossings":
        candidates = candidates.sort_values(["num_crossings", "id_num"])

    already_done = total_in_range - len(candidates)
    print(
        f"Starting {args.kind} worker: {len(candidates)} to compute, "
        f"{already_done} already done and skipped, "
        f"ids {args.start}-{end}, max_crossings={args.max_crossings}, "
        f"order={args.order}",
        flush=True,
    )

    if args.dry_run:
        preview = candidates[["id_num", "num_crossings"]].head(20)
        for id_num, crossings in preview.itertuples(index=False):
            print(f"  would compute id_num={int(id_num)} ({int(crossings)} crossings)")
        if len(candidates) > len(preview):
            print(f"  ... and {len(candidates) - len(preview)} more")
        return 0
    failures = 0
    for position, id_num in enumerate(candidates["id_num"].astype(int), start=1):
        print(f"[{position}/{len(candidates)}] id_num={id_num}", flush=True)
        try:
            if args.kind == "floer":
                compute_floer_homology(
                    id_num, args.max_crossings, data_path=path
                )
            else:
                compute_s_invariants(
                    id_num,
                    args.max_crossings,
                    args.timeout,
                    data_path=path,
                    jar_path=args.jar,
                )
        except Exception:
            failures += 1
            traceback.print_exc()
            print(f"Continuing after failure for id_num={id_num}", flush=True)

    print(f"Worker complete: failures={failures}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
