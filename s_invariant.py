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


KNOTJOB_JAR = (
    Path(__file__).resolve().parent / "tools" / "knotjob" / "KnotJob.jar"
)
_RESULT_PATTERN = re.compile(
    r"S-Invariant mod (\d+)\s*:\s*(-?\d+)", re.IGNORECASE
)


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


def _read_csv(
    csv_path: Path,
) -> tuple[list[str], list[str], list[dict[str, str]]]:
    """Read a CSV while retaining comments before its header."""
    lines = csv_path.read_text(encoding="utf-8").splitlines(keepends=True)
    header_index = next(
        (
            index
            for index, line in enumerate(lines)
            if line.strip() and not line.lstrip().startswith("#")
        ),
        None,
    )
    if header_index is None:
        raise ValueError(f"No CSV header was found in {csv_path}.")

    reader = csv.DictReader(io.StringIO("".join(lines[header_index:])))
    if reader.fieldnames is None:
        raise ValueError(f"No CSV header was found in {csv_path}.")
    rows = list(reader)
    if any(None in row for row in rows):
        raise ValueError("A CSV row has more values than the header.")
    return lines[:header_index], list(reader.fieldnames), rows


def _write_csv(
    csv_path: Path,
    preamble: Sequence[str],
    fieldnames: Sequence[str],
    rows: Sequence[dict[str, str]],
) -> None:
    """Atomically write the current CSV state."""
    file_mode = csv_path.stat().st_mode & 0o7777
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            dir=csv_path.parent,
            prefix=f".{csv_path.name}.",
            suffix=".tmp",
            encoding="utf-8",
            newline="",
            delete=False,
        ) as output:
            temporary_path = Path(output.name)
            output.writelines(preamble)
            if preamble and not preamble[-1].endswith(("\n", "\r")):
                output.write("\n")
            writer = csv.DictWriter(
                output,
                fieldnames=fieldnames,
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
            output.flush()
            os.fsync(output.fileno())
        os.chmod(temporary_path, file_mode)
        os.replace(temporary_path, csv_path)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def _s_invariant_column(prime: int) -> str:
    return f"s_invariant_mod_{prime}"


def _is_missing(value: str | None) -> bool:
    return value is None or value.strip().lower() in {"", "na", "n/a", "nan"}


def update_csv_s_invariants(
    csv_path: str | os.PathLike[str],
    primes: Iterable[int] = (2, 3),
    *,
    max_heap: str = "8g",
    java_path: str | os.PathLike[str] | None = None,
    jar_path: str | os.PathLike[str] = KNOTJOB_JAR,
    timeout: float | None = 10,
) -> None:
    """Calculate missing s-invariants and save each completed CSV row."""
    if timeout is not None and timeout <= 0:
        raise ValueError("timeout must be positive or None.")
    path = Path(csv_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"CSV file was not found: {path}")

    characteristics = _validate_primes(primes)
    preamble, fieldnames, rows = _read_csv(path)
    if "PD_code" not in fieldnames:
        raise ValueError("The CSV must contain a PD_code column.")

    columns = {
        prime: _s_invariant_column(prime) for prime in characteristics
    }
    changed = False
    for column in columns.values():
        if column not in fieldnames:
            fieldnames.append(column)
            changed = True
        for row in rows:
            if _is_missing(row.get(column)):
                if row.get(column) != "na":
                    changed = True
                row[column] = "na"

    if changed:
        _write_csv(path, preamble, fieldnames, rows)

    total = len(rows)
    for index, row in enumerate(rows, start=1):
        missing = [
            prime
            for prime, column in columns.items()
            if _is_missing(row[column])
        ]
        if not missing:
            continue

        try:
            pd_code = ast.literal_eval(row["PD_code"])
        except (SyntaxError, ValueError) as error:
            raise ValueError(f"Invalid PD_code in CSV row {index}.") from error

        knot_name = row.get("name") or "knot"
        label = f"{knot_name}_row_{index}"
        for prime in missing:
            print(
                f"[{index}/{total}] Calculating {label} mod {prime}...",
                flush=True,
            )
            try:
                results = calculate_s_invariants(
                    pd_code,
                    primes=(prime,),
                    name=label,
                    max_heap=max_heap,
                    java_path=java_path,
                    jar_path=jar_path,
                    timeout=timeout,
                )
            except KnotJobTimeoutError:
                print(
                    f"[{index}/{total}] Timed out mod {prime}; leaving na.",
                    flush=True,
                )
                continue

            value = results[prime]
            row[columns[prime]] = str(value)
            _write_csv(path, preamble, fieldnames, rows)
            print(
                f"[{index}/{total}] Saved mod {prime}: {value}.",
                flush=True,
            )


def _parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add missing Rasmussen s-invariants to a friends CSV."
    )
    parser.add_argument("csv_path", help="Path to a CSV containing PD_code.")
    parser.add_argument(
        "--primes",
        nargs="+",
        type=int,
        default=[2, 3],
        help="Coefficient characteristics to calculate (default: 2 3).",
    )
    parser.add_argument(
        "--max-heap",
        default="8g",
        help="Maximum Java heap size (default: 8g).",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=10,
        help="Timeout in seconds for each invariant (default: 10).",
    )
    return parser.parse_args()


def main() -> None:
    arguments = _parse_arguments()
    update_csv_s_invariants(
        arguments.csv_path,
        primes=arguments.primes,
        max_heap=arguments.max_heap,
        timeout=arguments.timeout,
    )


if __name__ == "__main__":
    main()
