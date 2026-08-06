"""Shared, concurrency-safe CSV helpers for invariant calculations."""

from __future__ import annotations

import fcntl
import os
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATA_PATH = PROJECT_ROOT / "Data" / "n_friends(master version).csv"


def resolve_data_path(path: str | os.PathLike[str] | None = None) -> Path:
    """Resolve the output CSV, allowing an environment-variable override."""
    configured = path or os.environ.get("N_FRIENDS_DATA") or DEFAULT_DATA_PATH
    return Path(configured).expanduser().resolve()


def value_is_missing(value: Any) -> bool:
    """Return whether a CSV invariant value has not been computed."""
    return pd.isna(value) or str(value).strip() in {"", "nan", "None"}


def update_row(
    id_num: int,
    updates: dict[str, Any],
    *,
    data_path: str | os.PathLike[str] | None = None,
) -> None:
    """Apply a row update without losing writes from another worker.

    The expensive invariant calculation happens before this function is called.
    Only the short read-modify-replace operation is protected by a file lock.
    """
    path = resolve_data_path(data_path)
    lock_path = path.with_suffix(path.suffix + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)

    with lock_path.open("a+") as lock_file:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        data = pd.read_csv(path)
        matches = data.index[data["id_num"] == int(id_num)]
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected one row with id_num={id_num}, found {len(matches)}."
            )
        row_index = matches[0]
        for column, value in updates.items():
            data.at[row_index, column] = value

        descriptor, temporary_name = tempfile.mkstemp(
            prefix=path.name + ".", suffix=".tmp", dir=path.parent
        )
        os.close(descriptor)
        temporary_path = Path(temporary_name)
        try:
            data.to_csv(temporary_path, index=False)
            os.replace(temporary_path, path)
        finally:
            temporary_path.unlink(missing_ok=True)
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)
