import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
TESTS = REPO / "tests"
FIX = TESTS / "fixtures"
FAKE_BIN = TESTS / "fake_bin"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(TESTS))

import make_fake_inputs  # noqa: E402


def small_config(cfg, monomers):
    """config restricted to some monomers, both templates and their complexes."""
    keep = set(monomers) | {"histamine", "histidine"}
    keep |= {f"{m}_{t}_1to1" for m in monomers for t in ("histamine", "histidine")}
    out = dict(cfg)
    out["monomers"] = [m for m in cfg["monomers"] if m in monomers]
    out["species"] = {k: v for k, v in cfg["species"].items() if k in keep}
    out["bindings"] = [b for b in cfg["bindings"] if b["complex"] in keep]
    out["selectivity"] = [s for s in cfg["selectivity"] if s[0].split()[0] in monomers]
    return out


class Run:
    """A scratch directory with pipeline.py, a config, fake inputs and fake programs."""

    def __init__(self, path, cfg):
        self.path = path
        shutil.copy(REPO / "pipeline.py", path / "pipeline.py")
        cfg = dict(cfg, crest="crest", orca=str(FAKE_BIN / "orca"))
        (path / "config.json").write_text(json.dumps(cfg, indent=1))
        make_fake_inputs.main(path / "config.json", path)
        self.env = {**os.environ,
                    "PATH": f"{FAKE_BIN}{os.pathsep}{os.environ['PATH']}",
                    "FAKE_SLURM_STATE": str(path / ".fake_slurm")}

    def __call__(self, *args, env=None, check=True):
        p = subprocess.run([sys.executable, "pipeline.py", *args],
                           cwd=self.path, env={**self.env, **(env or {})},
                           capture_output=True, text=True)
        if check and p.returncode != 0:
            raise AssertionError(f"pipeline {args} failed:\n{p.stdout}\n{p.stderr}")
        return p

    def tool(self, name):
        return subprocess.run([str(FAKE_BIN / name)], cwd=self.path, env=self.env,
                              capture_output=True, text=True, check=True)


@pytest.fixture
def water_cfg():
    return json.loads((REPO / "config_water.json").read_text())


@pytest.fixture
def make_run(tmp_path):
    def make(cfg):
        return Run(tmp_path, cfg)
    return make


def parse_status(text):
    """{(species, stage): {state: n}} from `pipeline.py status` output."""
    states = ["done", "running", "queued", "failed", "not started"]
    rows, species = {}, None
    for line in text.splitlines():
        parts = line.split()
        if len(parts) < 6 or not parts[-1].isdigit():
            continue
        if len(parts) == 7:
            species = parts[0]
            parts = parts[1:]
        rows[(species, parts[0])] = dict(zip(states, map(int, parts[1:])))
    return rows
