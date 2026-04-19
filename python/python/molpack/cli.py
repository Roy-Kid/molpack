"""molpack command-line interface.

The ``molpack pack`` subcommand reads a molpack ``.inp`` script — the
same line-oriented, keyword-driven format the Rust binary consumes —
and runs the packing loop.

Example::

    $ cat mixture.inp
    tolerance 2.0
    seed 1234567
    filetype pdb
    output mixture.pdb

    structure water.pdb
      number 1000
      inside box 0. 0. 0. 40. 40. 40.
    end structure

    structure urea.pdb
      number 400
      inside box 0. 0. 0. 40. 40. 40.
    end structure

    $ molpack pack mixture.inp
"""

from __future__ import annotations

import sys
from importlib.metadata import PackageNotFoundError, version as _pkg_version
from pathlib import Path
from typing import Annotated, Any

import typer

from . import PackError, load_script

app = typer.Typer(
    name="molpack",
    help="Molecular packing from a `.inp` script.",
    no_args_is_help=True,
    add_completion=False,
)


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def _require_molrs() -> Any:
    try:
        import molrs
    except ImportError as exc:
        raise typer.BadParameter(
            "molrs is required to read/write molecular files; install 'molcrafts-molrs'"
        ) from exc
    return molrs


def _read_frame(path: Path) -> Any:
    if not path.exists():
        raise typer.BadParameter(f"structure file not found: {path}")
    suffix = path.suffix.lower()
    molrs = _require_molrs()
    if suffix == ".pdb":
        return molrs.read_pdb(str(path))
    if suffix == ".xyz":
        return molrs.read_xyz(str(path))
    raise typer.BadParameter(f"unsupported input extension '{suffix}' for {path}")


def _write_result(path: Path, result: Any) -> None:
    molrs = _require_molrs()
    import numpy as np

    frame_dict = result.frame
    atoms = frame_dict["atoms"]
    block = molrs.Block()
    block.insert("element", list(atoms["element"]))
    block.insert("x", np.asarray(atoms["x"], dtype=np.float64))
    block.insert("y", np.asarray(atoms["y"], dtype=np.float64))
    block.insert("z", np.asarray(atoms["z"], dtype=np.float64))
    frame = molrs.Frame()
    frame["atoms"] = block

    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix.lower()
    if suffix == ".pdb":
        molrs.write_pdb(str(path), frame)
    elif suffix == ".xyz":
        molrs.write_xyz(str(path), frame)
    else:
        raise typer.BadParameter(f"unsupported output extension '{suffix}' for {path}")


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------


def _version_string() -> str:
    try:
        return _pkg_version("molcrafts-molpack")
    except PackageNotFoundError:
        return "unknown"


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"molpack {_version_string()}")
        raise typer.Exit()


@app.callback()
def _main(
    version: Annotated[
        bool,
        typer.Option(
            "--version",
            "-V",
            callback=_version_callback,
            is_eager=True,
            help="Show version and exit.",
        ),
    ] = False,
) -> None:
    """Packmol-grade molecular packing."""


@app.command()
def version() -> None:
    """Print the installed molpack version."""
    typer.echo(f"molpack {_version_string()}")


@app.command()
def info(
    path: Annotated[Path, typer.Argument(help="Path to a .pdb or .xyz file.")],
) -> None:
    """Print atom counts and bounding box for a structure file."""
    frame = _read_frame(path)
    atoms = frame["atoms"]
    try:
        elements = list(atoms.view("element"))
    except Exception:
        elements = list(atoms.view("symbol"))
    xs = list(atoms.view("x"))
    ys = list(atoms.view("y"))
    zs = list(atoms.view("z"))
    counts: dict[str, int] = {}
    for el in elements:
        counts[el] = counts.get(el, 0) + 1
    typer.echo(f"file     : {path}")
    typer.echo(f"natoms   : {len(elements)}")
    breakdown = ", ".join(f"{e}={n}" for e, n in sorted(counts.items()))
    typer.echo(f"elements : {breakdown}")
    if xs:
        typer.echo(
            f"extent   : x[{min(xs):.3f},{max(xs):.3f}] "
            f"y[{min(ys):.3f},{max(ys):.3f}] "
            f"z[{min(zs):.3f},{max(zs):.3f}]"
        )


@app.command()
def pack(
    script: Annotated[Path, typer.Argument(help="Path to a molpack `.inp` script.")],
    output: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help="Output file (.pdb or .xyz). Overrides the script's `output` keyword.",
        ),
    ] = None,
    seed: Annotated[
        int | None,
        typer.Option(
            "--seed", "-s", help="Random seed. Overrides the script's `seed` keyword."
        ),
    ] = None,
    max_loops: Annotated[
        int | None,
        typer.Option(
            "--max-loops",
            "-n",
            help="Maximum outer-loop iterations. Overrides the script's `nloop`.",
        ),
    ] = None,
    progress: Annotated[
        bool,
        typer.Option(
            "--progress/--no-progress",
            help="Show/suppress the built-in progress handler.",
        ),
    ] = True,
) -> None:
    """Run a packing job from a `.inp` script."""
    if not script.exists():
        typer.echo(f"error: script not found: {script}", err=True)
        raise typer.Exit(code=2)

    job = load_script(script)

    packer = job.packer.with_progress(progress)
    if seed is not None:
        packer = packer.with_seed(int(seed))

    effective_output = output if output is not None else Path(job.output)
    effective_max_loops = max_loops if max_loops is not None else job.nloop

    typer.echo(
        f"packing {len(job.targets)} target(s), "
        f"{sum(t.count for t in job.targets)} molecule(s), "
        f"max_loops={effective_max_loops}"
    )

    try:
        result = packer.pack(job.targets, max_loops=effective_max_loops)
    except PackError as exc:
        typer.echo(f"error: packing failed: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    typer.echo(
        f"done: converged={result.converged} natoms={result.natoms} "
        f"fdist={result.fdist:.4f} frest={result.frest:.4f}"
    )

    _write_result(effective_output, result)
    typer.echo(f"wrote {effective_output}")

    if not result.converged:
        raise typer.Exit(code=1)


def main() -> None:
    try:
        app()
    except KeyboardInterrupt:
        typer.echo("interrupted", err=True)
        sys.exit(130)


if __name__ == "__main__":
    main()
