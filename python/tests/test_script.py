import pathlib

import molrs
import numpy as np

import molpack

WATER_XYZ = """3
water
O 0.000 0.000 0.000
H 0.957 0.000 0.000
H -0.240 0.927 0.000
"""

SCRIPT = """tolerance 2.0
seed 7
output packed.xyz

structure water.xyz
  number 4
  inside cube 0. 0. 0. 12.
end structure
"""


def _write_inputs(tmp_path: pathlib.Path) -> pathlib.Path:
    (tmp_path / "water.xyz").write_text(WATER_XYZ)
    script = tmp_path / "mix.inp"
    script.write_text(SCRIPT)
    return script


def test_load_script_default_loader_reads_templates_and_packs(tmp_path):
    job = molpack.load_script(_write_inputs(tmp_path))
    packer, targets, output, nloop = job
    assert job.packer is packer
    assert [t.natoms for t in targets] == [3]
    assert [t.count for t in targets] == [4]
    assert output == tmp_path.resolve() / "packed.xyz"
    assert nloop == 200
    state = packer.run(targets, max_loops=nloop)
    assert state.frame["atoms"]["x"].shape == (12,)


def test_load_script_loader_replaces_the_molrs_reader(tmp_path):
    seen: list[tuple[str, str | None]] = []

    def loader(path: str, filetype: str | None) -> molrs.core.Frame:
        seen.append((path, filetype))
        xyz = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        return molrs.core.Frame(
            {
                "atoms": {
                    "x": xyz[:, 0].copy(),
                    "y": xyz[:, 1].copy(),
                    "z": xyz[:, 2].copy(),
                    "element": ["O", "H"],
                }
            }
        )

    job = molpack.load_script(_write_inputs(tmp_path), loader=loader)
    assert seen == [(str(tmp_path.resolve() / "water.xyz"), None)]
    assert [t.natoms for t in job.targets] == [2]
