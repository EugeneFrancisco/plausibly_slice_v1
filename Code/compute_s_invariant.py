"""Compute Rasmussen s-invariants from planar diagram codes."""

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
import caffeine
from tqdm import tqdm
import math


#Change these file paths to match your own computer
KNOTJOB_JAR = Path("/Users/henrigreamo/Desktop/KnotJob/KnotJob.jar")

out_file_path = "/Users/henrigreamo/Desktop/plausibly_slice_v1/Data/n_friends.csv"

load('/Users/henrigreamo/Desktop/plausibly_slice_v1/Code/n_friends_search.py')

_RESULT_PATTERN = re.compile(
    r"S-Invariant mod (\d+)\s*:\s*(-?\d+)", re.IGNORECASE
)

data_types = {
    "id_num": int,
    "num_crossings": int,
    "volume": float,
    "n": int,
    "verification": bool,
    "s": str,
    "s_2": str,
    "s_3": str,
    "obstructs": str,
    "n_friend_name": str,
    "n_friend_index": int,
    "knot_PD_code": str,
    "n_friend_PD_code": str
    }

class KnotJobTimeoutError(RuntimeError):
    """Raised when KnotJob exceeds its calculation timeout."""


def _find_java(java_path: str | os.PathLike[str] | None = None) -> str:
    """Find the Java executable used to run KnotJob."""
    if java_path is not None:
        return os.fspath(java_path)

    configured = os.environ.get("KNOTJOB_JAVA")
    if configured:
        return configured

    for conda_root in ("miniforge3", "miniconda3", "anaconda3"):
        environment = Path.home() / conda_root / "envs" / "jdk23"
        for relative_path in ("bin/java", "lib/jvm/bin/java"):
            candidate = environment / relative_path
            if candidate.is_file():
                return str(candidate)

    java = shutil.which("java")
    if java:
        return java
    raise RuntimeError("Java was not found; set KNOTJOB_JAVA to Java 23 or newer.")


def _validate_primes(primes: Iterable[int]) -> tuple[int, ...]:
    """Validate and deduplicate KnotJob coefficient characteristics."""
    result = []
    for prime in primes:
        if isinstance(prime, bool) or not isinstance(prime, int):
            raise TypeError("Each characteristic must be an integer.")
        if prime != 0 and (
            prime >= 212
            or prime < 2
            or any(
                prime % divisor == 0
                for divisor in range(2, int(prime**0.5) + 1)
            )
        ):
            raise ValueError(
                "Each characteristic must be 0 or a prime less than 212."
            )
        if prime not in result:
            result.append(prime)
    if not result:
        raise ValueError("At least one characteristic is required.")
    return tuple(result)


def pd_code_to_knotjob(
    pd_code: Sequence[Sequence[int]], name: str = "knot"
) -> str:
    """Convert a PD code into KnotJob's text input format."""
    if not name or "\n" in name or "\r" in name or "=" in name:
        raise ValueError(
            "The knot name must be nonempty and contain no newline or '='."
        )

    crossings = []
    for crossing in pd_code:
        if len(crossing) != 4:
            raise ValueError("Each PD crossing must contain exactly four labels.")
        if any(
            isinstance(label, bool) or not isinstance(label, int)
            for label in crossing
        ):
            raise TypeError("PD labels must be integers.")
        crossings.append(f"X[{','.join(str(label) for label in crossing)}]")

    if not crossings:
        raise ValueError("The PD code must contain at least one crossing.")
    return f"{name} = PD {', '.join(crossings)}"


def calculate_s_invariants(
    pd_code: Sequence[Sequence[int]],
    primes: Iterable[int] = (2, 3),
    *,
    name: str = "knot",
    max_heap: str = "8g",
    java_path: str | os.PathLike[str] | None = None,
    jar_path: str | os.PathLike[str] = KNOTJOB_JAR,
    timeout: float | None = None,
) -> dict[int, int]:
    """Return Rasmussen s-invariants keyed by coefficient characteristic."""
    if timeout is not None and timeout <= 0:
        raise ValueError("timeout must be positive or None.")
    characteristics = _validate_primes(primes)
    knotjob_input = pd_code_to_knotjob(pd_code, name)
    jar = Path(jar_path)
    if not jar.is_file():
        raise RuntimeError(f"KnotJob.jar was not found at {jar}.")
    if not re.fullmatch(r"\d+[kKmMgG]?", max_heap):
        raise ValueError(
            "max_heap must be a Java heap size such as '8g' or '1024m'."
        )

    command = [
        _find_java(java_path),
        f"-Xmx{max_heap}",
        "-jar",
        str(jar),
    ]

    with tempfile.TemporaryDirectory() as temporary_directory:
        input_path = Path(temporary_directory) / "knot.txt"
        input_path.write_text(knotjob_input + "\n", encoding="utf-8")
        command.append(str(input_path))
        command.extend(f"-s{prime}" for prime in characteristics)
        command.append("-nf")
        try:
            process = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
        except FileNotFoundError as error:
            raise RuntimeError(
                f"Java executable was not found: {command[0]}"
            ) from error
        except subprocess.TimeoutExpired as error:
            raise KnotJobTimeoutError(
                "KnotJob exceeded the requested timeout."
            ) from error

    if process.returncode != 0:
        details = (process.stderr or process.stdout).strip()
        raise RuntimeError(
            f"KnotJob failed with exit code {process.returncode}: {details}"
        )

    results = {
        int(prime): int(value)
        for prime, value in _RESULT_PATTERN.findall(process.stdout)
    }
    missing = [prime for prime in characteristics if prime not in results]
    if missing:
        raise RuntimeError(
            f"KnotJob returned no s-invariant for characteristics {missing}. "
            f"Output: {process.stdout.strip()}"
        )
    return {prime: results[prime] for prime in characteristics}

def write_s_invariants(id_num: int, max_crossings: int = 45):
    data = pd.read_csv(out_file_path,dtype=data_types)
    row = data.iloc[int(id_num) - 1]

    #print(row)
    num_crossings = int(row["num_crossings"])
    pd_code = ast.literal_eval(row["knot_PD_code"])
    n = int(row["n"])

    if num_crossings > max_crossings:
        print("Too many crossings")
        return

    print(f"Finding s-invariants for knot {id_num}")
    results = calculate_s_invariants(pd_code=pd_code,primes=(0,2,3))
    print(results)
    row["s"]=str(results[0])
    row["s_2"]=str(results[2])
    row["s_3"]=str(results[3])

    obstructs = False
    for s in results.values():
        if (s > n-math.sqrt(n)):
            obstructs = True

    row["obstructs"]=str(obstructs)

    data.iloc[int(id_num)-1]=row
    data.to_csv(out_file_path,index=False)
    
def compute_specified_knots(indices: list[int]):
    for i in indices:
        write_s_invariants(i)


