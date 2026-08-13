"""Simplify knot_PD_code entries in the n-friends master CSV in place.

Standalone, plain-Python (no Sage). For each row it runs snappy's global
simplifier on knot_PD_code, verifies the simplified diagram has the same
hyperbolic volume, and writes the shorter PD code + updated num_crossings back
through invariant_io.update_row so the write stays atomic under the same lock
the invariant workers use.

A `simplified` column tracks which rows have already been processed so reruns
skip them. Progress is persisted after every row, so an interrupted run loses
at most one knot's worth of work.
"""

from __future__ import annotations

import argparse
import ast
import fcntl
import sys
import traceback
from pathlib import Path

import pandas as pd
import snappy
from tqdm import tqdm

from invariant_io import resolve_data_path, update_row


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default=None, help="n-friends CSV (defaults to master)")
    parser.add_argument("--start", type=int, default=1, help="first id_num to inspect")
    parser.add_argument("--end", type=int, default=None, help="last id_num to inspect")
    parser.add_argument(
        "--n",
        type=int,
        nargs="+",
        default=[4, 5],
        help="only process rows whose `n` column is in this list (default: 4 5)",
    )
    parser.add_argument(
        "--max-crossings",
        type=int,
        default=1000,
        help="skip rows whose current num_crossings exceeds this",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=5,
        help="outer iterations of super_simplify",
    )
    parser.add_argument("--type-3-limit", type=int, default=5000)
    parser.add_argument(
        "--volume-tol",
        type=float,
        default=1e-6,
        help="acceptable |new_volume - stored_volume| before we refuse to write",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="report what would be simplified, then exit without computing",
    )
    return parser.parse_args()


def ensure_simplified_column(path: Path) -> None:
    """Add the `simplified` column with default False if it doesn't exist."""
    lock_path = path.with_suffix(path.suffix + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as lock_file:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        try:
            data = pd.read_csv(path)
            if "simplified" not in data.columns:
                data["simplified"] = False
                tmp = path.with_suffix(path.suffix + ".bootstrap.tmp")
                data.to_csv(tmp, index=False)
                tmp.replace(path)
                print("Added `simplified` column (default False) to CSV", flush=True)
        finally:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


def super_simplify(link, iterations: int, type_3_limit: int):
    """Loop K.simplify('global', ...) until nothing reduces, N outer times."""
    K = link.copy()
    for _ in range(iterations):
        reduced = True
        while reduced:
            reduced = K.simplify("global", type_3_limit)
    return K


def already_simplified(value) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    return text in {"true", "1", "1.0", "t", "yes"}


def main() -> int:
    args = parse_args()
    path = resolve_data_path(args.data)
    ensure_simplified_column(path)
    data = pd.read_csv(path)

    end = args.end if args.end is not None else int(data["id_num"].max())
    n_filter = set(args.n)
    in_range = data.loc[data["id_num"].between(args.start, end)]

    to_process = []
    skipped_already = 0
    skipped_max = 0
    skipped_n = 0
    for _, row in in_range.iterrows():
        if int(row["n"]) not in n_filter:
            skipped_n += 1
            continue
        if already_simplified(row.get("simplified")):
            skipped_already += 1
            continue
        if int(row["num_crossings"]) > args.max_crossings:
            skipped_max += 1
            continue
        to_process.append(row)

    print(
        f"Simplify worker: {len(to_process)} to process, "
        f"{skipped_already} already simplified, "
        f"{skipped_max} above --max-crossings={args.max_crossings}, "
        f"{skipped_n} with n not in {sorted(n_filter)}, "
        f"ids {args.start}-{end}",
        flush=True,
    )

    if args.dry_run:
        for row in to_process[:20]:
            print(
                f"  would simplify id_num={int(row['id_num'])} "
                f"({int(row['num_crossings'])} crossings)"
            )
        if len(to_process) > 20:
            print(f"  ... and {len(to_process) - 20} more")
        return 0

    failures = 0
    reduced_count = 0
    unchanged_count = 0
    volume_mismatches = 0
    progress = tqdm(to_process, desc="Simplifying", unit="knot")
    for row in progress:
        id_num = int(row["id_num"])
        old_crossings = int(row["num_crossings"])
        expected_volume = float(row["volume"])
        try:
            pd_code = ast.literal_eval(row["knot_PD_code"])
            K = snappy.Link(pd_code)
            K = super_simplify(K, args.iterations, args.type_3_limit)
            new_pd = K.PD_code()
            new_crossings = len(K.crossings)

            try:
                new_volume = float(snappy.Link(new_pd).exterior().volume())
            except Exception:
                new_volume = None

            if new_volume is None or abs(new_volume - expected_volume) > args.volume_tol:
                volume_mismatches += 1
                tqdm.write(
                    f"id_num={id_num}: volume mismatch "
                    f"(expected {expected_volume}, got {new_volume}); skipping write"
                )
                progress.set_postfix(reduced=reduced_count, same=unchanged_count,
                                     vol_mismatch=volume_mismatches, fail=failures)
                continue

            updates = {"simplified": True}
            if new_crossings < old_crossings:
                updates["knot_PD_code"] = str(new_pd)
                updates["num_crossings"] = new_crossings
                reduced_count += 1
                tqdm.write(f"id_num={id_num}: reduced {old_crossings} -> {new_crossings}")
            else:
                unchanged_count += 1

            update_row(id_num, updates, data_path=path)
        except Exception:
            failures += 1
            tqdm.write(f"id_num={id_num}: FAILED")
            traceback.print_exc()
        progress.set_postfix(reduced=reduced_count, same=unchanged_count,
                             vol_mismatch=volume_mismatches, fail=failures)
    progress.close()

    print(
        f"Done: reduced={reduced_count}, unchanged={unchanged_count}, "
        f"volume_mismatches={volume_mismatches}, failures={failures}",
        flush=True,
    )
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
