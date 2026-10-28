"""pipeline.py end to end against the fake crest and orca in tests/fake_bin."""
import csv
import json
import math
import os
import shutil

import pytest

from conftest import FIX, parse_status, small_config

RT = 0.0019872043 * 298.15


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.reader(f))


def test_all_water_end_to_end(make_run, water_cfg):
    run = make_run(water_cfg)
    p = run("all", "config.json")
    res = run.path / "results"

    summary = read_csv(res / "summary.csv")
    assert summary[0] == ["level", "type", "name", "dG_kcal_1M", "K_ratio"]
    rows = summary[1:]
    nb, ns = len(water_cfg["bindings"]), len(water_cfg["selectivity"])
    assert (nb, ns) == (16, 8)
    assert len(rows) == 2 * (nb + ns)
    for level in ("screen", "high"):
        mine = [r for r in rows if r[0] == level]
        assert [r[2] for r in mine if r[1] == "binding"] == [b["name"] for b in water_cfg["bindings"]]
        gaps = [r for r in mine if r[1] == "selectivity gap (C1)"]
        assert len(gaps) == ns
        for r in mine:
            assert -50 < float(r[3]) < 50
        for r in gaps:
            assert float(r[4]) == pytest.approx(math.exp(-float(r[3]) / RT), rel=0.02, abs=0.01)

    ranking = read_csv(res / "ranking.csv")
    assert ranking[0] == ["monomer", "dG_bind_histamine_kcal", "selectivity_gap_kcal",
                          "oxidation_potential_V_vs_SHE", "measured_IF"]
    assert [r[0] for r in ranking[1:]] == water_cfg["monomers"]
    high = {r[2]: float(r[3]) for r in rows if r[0] == "high" and r[1] == "binding"}
    for m, dg, gap, ox, measured in ranking[1:]:
        assert float(dg) == pytest.approx(high[f"{m} histamine"], abs=0.006)
        assert float(gap) == pytest.approx(high[f"{m} histamine"] - high[f"{m} histidine"], abs=0.011)
        assert 0 < float(ox) < 10
        assert measured == ""

    for level in ("screen", "high"):
        confs = read_csv(res / f"conformers_{level}.csv")
        assert confs[0] == ["species", "variant", "conf", "dG_rel_kcal", "weight"]
        per = {}
        for sp, var, conf, dg, w in confs[1:]:
            per.setdefault(sp, []).append(float(w))
            assert float(dg) >= 0
        assert set(per) == set(water_cfg["species"])
        for sp, ws in per.items():
            assert sum(ws) == pytest.approx(1, abs=0.01)

    st = parse_status(run("status", "config.json").stdout)
    for stage in ("crest", "sp", "opt", "hl", "extras"):
        tot = st[("TOTAL", stage)]
        assert tot["done"] > 0
        assert tot["running"] + tot["queued"] + tot["failed"] + tot["not started"] == 0
    assert st[("TOTAL", "crest")]["done"] == sum(len(s["variants"]) for s in water_cfg["species"].values())
    assert "excluded (imaginary frequency)" in p.stdout

    # every stage is resumable: a second pass starts nothing
    again = run("all", "config.json")
    assert "  ORCA " not in again.stdout
    assert "[crest] " in again.stdout and "(done)" in again.stdout
    assert read_csv(res / "summary.csv") == summary


def test_status_counts_failures_and_progress(make_run, water_cfg):
    run = make_run(small_config(water_cfg, ["MAA", "HEMA"]))
    env = {"FAKE_CREST_FAIL": "HEMA_histidine_1to1/ionpair",
           "FAKE_ORCA_FAIL": "MAA/neutral/sp/c001",
           "FAKE_ORCA_KILL": "/histamine/cation/opt/"}
    run("crest", "config.json", env=env)
    run("sp", "config.json", env=env)

    st = parse_status(run("status", "config.json").stdout)
    assert st[("HEMA_histidine_1to1", "crest")]["failed"] == 1
    assert st[("HEMA_histidine_1to1", "crest")]["done"] == 2
    assert st[("MAA", "sp")]["failed"] == 1
    assert st[("MAA", "opt")]["not started"] > 0
    assert ("MAA", "hl") not in st                  # nothing known yet

    out = run("opt", "config.json", env=env).stdout
    assert "!! failed" in out
    text = run("status", "config.json").stdout
    st = parse_status(text)
    assert st[("histamine", "opt")]["failed"] >= 1
    # every histamine opt is a cation (see REVIEW.md, pooled charge states), so no hl
    assert st[("histamine", "opt")]["done"] == 0
    assert ("histamine", "hl") not in st
    assert st[("MAA", "hl")]["not started"] >= 1
    assert "failed jobs (" in text
    assert "work/MAA/neutral/sp/c001/sp.out" in text
    assert "work/HEMA_histidine_1to1/ionpair/crest/crest.log" in text

    # a rerun without the failure switches heals everything
    run("crest", "config.json")
    run("sp", "config.json")
    run("opt", "config.json")
    st = parse_status(run("status", "config.json").stdout)
    for stage in ("crest", "sp", "opt"):
        assert st[("TOTAL", stage)]["failed"] == 0


def test_running_marker(make_run, water_cfg):
    run = make_run(small_config(water_cfg, ["MAA"]))
    run("crest", "config.json")
    run("sp", "config.json")
    d = run.path / "work" / "MAA" / "neutral" / "sp" / "c000"
    (d / "sp.out").write_text("half an output\n")
    host = os.uname().nodename

    # a live process holds the job: status says running and sp leaves it alone
    (d / "sp.running").write_text(f"local {os.getpid()} {host}\n")
    st = parse_status(run("status", "config.json").stdout)
    assert st[("MAA", "sp")]["running"] == 1
    out = run("sp", "config.json").stdout
    assert "(already running)" in out
    assert (d / "sp.out").read_text() == "half an output\n"
    assert "running or queued, skipping" in run("opt", "config.json").stdout

    # the process is gone: the job counts as failed and a rerun redoes it
    (d / "sp.running").write_text(f"local 999999999 {host}\n")
    st = parse_status(run("status", "config.json").stdout)
    assert st[("MAA", "sp")]["failed"] == 1
    run("sp", "config.json")
    assert "TERMINATED NORMALLY" in (d / "sp.out").read_text()
    assert not (d / "sp.running").exists()


def test_interrupted_crest_is_not_done(make_run, water_cfg):
    run = make_run(small_config(water_cfg, ["MAA"]))
    run("crest", "config.json")
    wd = run.path / "work" / "MAA" / "neutral" / "crest"
    mid = FIX / "crest3" / "maa_water_midrun"
    shutil.copy(mid / "crest.log", wd / "crest.log")
    shutil.copy(mid / "crest_conformers.xyz", wd / "crest_conformers.xyz")

    st = parse_status(run("status", "config.json").stdout)
    assert st[("MAA", "crest")]["failed"] == 1
    assert "no finished CREST output for MAA/neutral" in run("sp", "config.json").stdout
    out = run("crest", "config.json").stdout
    assert "[crest] MAA/neutral\n" in out
    assert parse_status(run("status", "config.json").stdout)[("MAA", "crest")]["done"] == 1


def test_missing_orca_is_a_clear_error(make_run, water_cfg):
    run = make_run(small_config(water_cfg, ["MAA"]))
    run("crest", "config.json")
    c = json.loads((run.path / "config.json").read_text())
    c["orca"] = "/nonexistent/orca"
    (run.path / "config.json").write_text(json.dumps(c))
    p = run("sp", "config.json", check=False)
    assert p.returncode != 0
    assert "cannot run /nonexistent/orca" in p.stderr


def test_extras_keeps_hand_labelled_cp_geometry(make_run, water_cfg):
    run = make_run(small_config(water_cfg, ["MAA"]))
    run("all", "config.json", env={"FAKE_ORCA_FAIL": "/extras/cp"})
    (geom,) = (run.path / "work" / "MAA_histamine_1to1").glob("*/extras/cp/geom.xyz")
    labelled = geom.read_text().replace("\n C ", "\n C(1) ", 1)
    assert "C(1)" in labelled
    geom.write_text(labelled)
    run("extras", "config.json")
    assert geom.read_text() == labelled
    assert "TERMINATED NORMALLY" in geom.with_name("cp.out").read_text()
