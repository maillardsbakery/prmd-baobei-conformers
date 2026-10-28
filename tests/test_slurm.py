"""--slurm against the fake sbatch/squeue in tests/fake_bin.

Fake sbatch records each array; its tasks stay pending (fake squeue lists
them) until fake-slurm-drain runs them, which is the scheduler's part.
"""
import csv
import json
import os
import re

from conftest import parse_status, small_config


def sbatch_files(run, stage):
    return sorted((run.path / "slurm").glob(f"{stage}-*.sbatch"))


def test_all_with_slurm_is_refused(make_run, water_cfg):
    run = make_run(small_config(water_cfg, ["MAA"]))
    p = run("all", "config.json", "--slurm", check=False)
    assert p.returncode != 0
    assert "one stage at a time" in p.stderr


def test_sp_job_array_lifecycle(make_run, water_cfg):
    run = make_run(small_config(water_cfg, ["MAA"]))
    run("crest", "config.json")
    p = run("sp", "config.json", "--slurm")
    assert "[slurm] submitted" in p.stdout
    (script,) = sbatch_files(run, "sp")
    text = script.read_text()
    jobs = (script.with_suffix(".jobs")).read_text().split()
    n = len(jobs)
    assert n > 10
    assert f"#SBATCH --array=1-{n}%50" in text
    assert "#SBATCH --ntasks=8\n" in text and "#SBATCH --cpus-per-task=1\n" in text
    assert "#SBATCH --mem-per-cpu=3750M" in text
    assert "#SBATCH --time=04:00:00" in text
    assert "--partition" not in text
    assert script.with_suffix(".jobid").read_text().strip() == "1000"
    for j in jobs:
        assert os.access(j, os.X_OK)
        assert j.endswith("/sp.sh")
        assert (run.path / j).parent.joinpath("sp.queued").exists()
        assert (run.path / j).parent.joinpath("sp.inp").exists()
        assert not (run.path / j).parent.joinpath("sp.out").exists()

    st = parse_status(run("status", "config.json").stdout)
    assert st[("TOTAL", "sp")]["queued"] == n

    # rerunning the stage does not submit the same jobs twice, and opt waits
    p = run("sp", "config.json", "--slurm")
    assert "(already queued)" in p.stdout
    assert "[slurm] submitted" not in p.stdout
    assert len(sbatch_files(run, "sp")) == 1
    p = run("opt", "config.json", "--slurm")
    assert "still running or queued, skipping" in p.stdout
    assert not sbatch_files(run, "opt")

    run.tool("fake-slurm-drain")
    st = parse_status(run("status", "config.json").stdout)
    assert st[("TOTAL", "sp")]["done"] == n
    assert not list((run.path / "work").rglob("*.queued"))
    assert not list((run.path / "work").rglob("*.running"))


def test_slurm_and_serial_give_same_results(make_run, water_cfg, tmp_path_factory):
    cfg = small_config(water_cfg, ["MAA", "AMPS"])
    serial = make_run(cfg)
    serial("all", "config.json")

    from conftest import Run
    other = tmp_path_factory.mktemp("slurm_run")
    sl = Run(other, cfg)
    for stage in ("crest", "sp", "opt", "hl", "extras"):
        p = sl(stage, "config.json", "--slurm")
        assert "[slurm] submitted" in p.stdout, (stage, p.stdout)
        sl.tool("fake-slurm-drain")
    sl("analyze", "config.json")

    crest_script = sbatch_files(sl, "crest")[0].read_text()
    assert "#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=8" in crest_script
    assert "export OMP_NUM_THREADS=8" in crest_script
    crest_job = next((other / "work").rglob("crest.sh")).read_text()
    assert 'export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"' in crest_job
    assert "--alpb water" in crest_job

    for f in ("summary.csv", "ranking.csv", "conformers_high.csv"):
        a = list(csv.reader(open(serial.path / "results" / f)))
        b = list(csv.reader(open(other / "results" / f)))
        assert a == b, f
    st = parse_status(sl("status", "config.json").stdout)
    for stage in ("crest", "sp", "opt", "hl", "extras"):
        assert st[("TOTAL", stage)]["done"] > 0
        assert st[("TOTAL", stage)]["failed"] == 0


def test_slurm_config_block_and_chunking(make_run, water_cfg):
    cfg = small_config(water_cfg, ["MAA"])
    cfg["slurm"] = {"partition": "chem", "account": "thesis42", "max_array_size": 4,
                    "max_parallel": 2, "time": {"sp": "00:30:00"},
                    "sbatch_lines": ["#SBATCH --qos=short"],
                    "setup": ["module load orca/6.0.1"]}
    run = make_run(cfg)
    run("crest", "config.json")
    run("sp", "config.json", "--slurm")
    scripts = sbatch_files(run, "sp")
    n_jobs = sum(len(s.with_suffix(".jobs").read_text().split()) for s in scripts)
    assert len(scripts) == -(-n_jobs // 4)
    for s in scripts:
        t = s.read_text()
        size = len(s.with_suffix(".jobs").read_text().split())
        assert size <= 4
        assert f"#SBATCH --array=1-{size}%2" in t
        for line in ("#SBATCH --partition=chem", "#SBATCH --account=thesis42",
                     "#SBATCH --time=00:30:00", "#SBATCH --qos=short",
                     "module load orca/6.0.1"):
            assert line in t
    ids = {s.with_suffix(".jobid").read_text().strip() for s in scripts}
    assert len(ids) == len(scripts)
    assert parse_status(run("status", "config.json").stdout)[("TOTAL", "sp")]["queued"] == n_jobs
    run.tool("fake-slurm-drain")
    assert parse_status(run("status", "config.json").stdout)[("TOTAL", "sp")]["done"] == n_jobs


def test_no_submit_and_missing_sbatch(make_run, water_cfg):
    cfg = small_config(water_cfg, ["MAA"])
    run = make_run(cfg)
    run("crest", "config.json")
    p = run("sp", "config.json", "--slurm", "--no-submit")
    assert "not submitted" in p.stdout
    assert sbatch_files(run, "sp")
    assert not list((run.path / "work").rglob("*.queued"))
    st = parse_status(run("status", "config.json").stdout)
    assert st[("TOTAL", "sp")]["not started"] > 0

    c = json.loads((run.path / "config.json").read_text())
    c["slurm"] = {"sbatch": "/nonexistent/sbatch"}
    (run.path / "config.json").write_text(json.dumps(c))
    p = run("sp", "config.json", "--slurm")
    assert "sbatch failed" in p.stdout and "submit it by hand" in p.stdout
    assert not list((run.path / "work").rglob("*.queued"))


def test_task_killed_by_slurm_counts_as_failed(make_run, water_cfg):
    run = make_run(small_config(water_cfg, ["MAA"]))
    run("crest", "config.json")
    run("sp", "config.json")
    d = run.path / "work" / "MAA" / "neutral" / "sp" / "c000"
    (d / "sp.out").write_text("cut off by the wall-time limit\n")
    # the job script never got to remove its marker; squeue no longer knows the task
    (d / "sp.running").write_text("slurm 4242_7 node17\n")
    text = run("status", "config.json").stdout
    st = parse_status(text)
    assert st[("MAA", "sp")]["failed"] == 1
    assert re.search(r"work/MAA/neutral/sp/c000/sp.out", text)
    p = run("sp", "config.json", "--slurm")
    assert "[slurm] submitted 1 jobs" in p.stdout
