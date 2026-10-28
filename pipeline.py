#!/usr/bin/env python3
"""
Baobei 2.0 conformer pipeline

  crest    : CREST (GFN2-xTB, ALPB solvent) conformer search per species/variant
  sp       : r2SCAN-3c single points on every CREST conformer inside sp_window
  opt      : pool all variants of a species, rank by SP energy, Opt+Freq the best
  hl       : high-level single points (second functional) on optimized geometries
  extras   : C4 descriptors: counterpoise interaction energy and oxidation potential
  analyze  : ensemble free energies, binding, selectivity gap (C1), C2/C3 ranking table
  all      : run every stage in order
  status   : per species and stage, how many jobs are done, running, queued,
             failed or not started

Hypotheses this tests (see "Thesis v2 hypotheses"):
  C1  selectivity gap > 0 for at least one monomer, in implicit water
  C2  binding free energy predicts sensor performance in rank order
  C3  or the relationship is a volcano: a middle binding strength wins
  C4  other DFT descriptors beat binding energy alone

C2 and C3 need several monomers and sensor data; this script produces the
computational column (results/ranking.csv) that you later correlate against
measured imprinting factors.

Every stage is resumable: finished jobs (ORCA TERMINATED NORMALLY / CREST
terminated normally) are skipped, so you can stop, add compute, and rerun.

With --slurm, a stage writes its jobs as SLURM job arrays under slurm/ and
submits them with sbatch instead of running them one after another. Run one
stage at a time, check `status`, and start the next stage once it shows the
previous one done. An optional "slurm" block in the config sets partition,
account, wall times and module loads (see SLURM_DEFAULTS).

Usage:  python pipeline.py <stage> [config.json] [--slurm [--no-submit]]
"""
import argparse
import csv
import getpass
import json
import math
import os
import re
import shlex
import shutil
import socket
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

H2KCAL = 627.509474          # Eh -> kcal/mol
R_KCAL = 0.0019872043        # kcal/(mol K)
WORK = Path("work")
RESULTS = Path("results")
SLURM_DIR = Path("slurm")

# Overridable with a "slurm" block in the config. Wall times are per job.
SLURM_DEFAULTS = {
    "sbatch": "sbatch",
    "squeue": "squeue",
    "partition": None,
    "account": None,
    "time": {"crest": "24:00:00", "sp": "04:00:00", "opt": "48:00:00",
             "hl": "12:00:00", "extras": "24:00:00"},
    "mem_per_cpu_mb": None,      # default: 1.25 x maxcore, so ORCA has headroom
    "max_parallel": 50,          # array throttle (%N)
    "max_array_size": 1000,      # split bigger stages; SLURM's MaxArraySize is often 1001
    "sbatch_lines": [],          # extra lines, e.g. "#SBATCH --qos=normal"
    "setup": [],                 # shell lines before each job, e.g. "module load orca/6.0.1"
}


# ---------------------------------------------------------------- helpers
def load_cfg(path):
    with open(path) as f:
        return json.load(f)


def variants(cfg):
    for sp, spec in cfg["species"].items():
        for var, v in spec["variants"].items():
            yield sp, spec, var, v


def read_multixyz(path):
    """Return [(energy_Eh, xyz_block)] from a CREST multi-structure xyz file."""
    lines = Path(path).read_text().splitlines()
    out, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        comment = lines[i + 1]
        atoms = lines[i + 2:i + 2 + n]
        m = re.search(r"-?\d+\.\d+", comment)
        e = float(m.group()) if m else None
        out.append((e, f"{n}\n{comment}\n" + "\n".join(atoms) + "\n"))
        i += 2 + n
    return out


def orca_done(out):
    out = Path(out)
    return out.exists() and "ORCA TERMINATED NORMALLY" in out.read_text(errors="ignore")


def final_energy(out):
    hits = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)",
                      Path(out).read_text(errors="ignore"))
    return float(hits[-1]) if hits else None


def g_minus_eel(out):
    m = re.findall(r"G-E\(el\)\s+\.\.\.\s+(-?\d+\.\d+)",
                   Path(out).read_text(errors="ignore"))
    return float(m[-1]) if m else None


def has_imaginary(out):
    return "***imaginary mode***" in Path(out).read_text(errors="ignore")


def crest_done(wd):
    """CREST 3 rewrites crest_conformers.xyz after every cycle, so the file alone
    does not mean the search finished: also require the termination line."""
    wd = Path(wd)
    log = wd / "crest.log"
    if not (wd / "crest_conformers.xyz").exists() or not log.exists():
        return False
    with open(log, "rb") as f:
        f.seek(max(0, log.stat().st_size - 65536))
        return b"CREST terminated normally" in f.read()


# ---------------------------------------------------------------- job state
# A job is one program run in one directory, named after its log (crest.log,
# sp.out, ...). While it runs, <name>.running says who runs it; a job handed to
# SLURM that has not started yet has <name>.queued. `status` and the stages
# read these so a rerun never starts a second copy of a live job.
_squeue = {}


def squeue_ids(cfg):
    """Set of 'jobid_task' strings SLURM knows about, or None if squeue is unusable."""
    if "ids" not in _squeue:
        cmd = [slurm_cfg(cfg)["squeue"], "-h", "-r", "-o", "%i", "-u", getpass.getuser()]
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            _squeue["ids"] = set(p.stdout.split()) if p.returncode == 0 else None
        except (OSError, subprocess.TimeoutExpired):
            _squeue["ids"] = None
    return _squeue["ids"]


def marker_alive(cfg, marker):
    words = marker.read_text().split()
    if not words:
        return False
    if words[0] == "local":                    # local <pid> <host>
        if words[2] != socket.gethostname():
            return True                        # can't see another machine's processes
        try:
            os.kill(int(words[1]), 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            pass
        return True
    if words[0] == "slurm":                    # slurm <jobid>_<task> <host>
        jid = words[1]
    elif words[0] == "queued":                 # queued <manifest> <task>
        jobid = Path(words[1]).with_suffix(".jobid")
        if not jobid.exists():
            return False
        jid = f"{jobid.read_text().strip()}_{words[2]}"
    else:
        return False
    ids = squeue_ids(cfg)
    return True if ids is None else jid in ids


def job_state(cfg, wd, name, done, log):
    """One of: done, running, queued, failed, not started."""
    wd = Path(wd)
    if done:
        return "done"
    for state in ("running", "queued"):
        m = wd / f"{name}.{state}"
        if m.exists() and marker_alive(cfg, m):
            return state
    return "failed" if (wd / log).exists() else "not started"


def orca_state(cfg, wd, name):
    return job_state(cfg, wd, name, orca_done(Path(wd) / f"{name}.out"), f"{name}.out")


def crest_state(cfg, wd):
    return job_state(cfg, wd, "crest", crest_done(wd), "crest.log")


def active(cfg, wd, name):
    """True if this job is running or waiting in the SLURM queue."""
    wd = Path(wd)
    return any((wd / f"{name}.{s}").exists() and marker_alive(cfg, wd / f"{name}.{s}")
               for s in ("running", "queued"))


# ---------------------------------------------------------------- job runners
def run(cmd, cwd, logname, env=None):
    name = Path(logname).stem
    marker = Path(cwd) / f"{name}.running"
    marker.write_text(f"local {os.getpid()} {socket.gethostname()}\n")
    full_env = None
    if env:
        full_env = dict(os.environ)
        for k, val in env.items():
            full_env.setdefault(k, val)
    try:
        with open(Path(cwd) / logname, "w") as log:
            return subprocess.run(cmd, cwd=cwd, stdout=log,
                                  stderr=subprocess.STDOUT, env=full_env).returncode
    finally:
        if marker.exists():
            marker.unlink()


def slurm_cfg(cfg):
    sc = {**SLURM_DEFAULTS, **cfg.get("slurm", {})}
    sc["time"] = {**SLURM_DEFAULTS["time"], **cfg.get("slurm", {}).get("time", {})}
    return sc


class Runner:
    """Runs jobs as they come, or with slurm=True collects them into job arrays."""

    def __init__(self, cfg=None, slurm=False, submit=True):
        self.cfg, self.slurm, self.submit = cfg, slurm, submit
        self.pending = []

    def add(self, wd, cmd, logname, env=None):
        """Serial: run now and return the exit code. SLURM: queue it, return None."""
        if not self.slurm:
            return run(cmd, wd, logname, env)
        self.pending.append((Path(wd).resolve(), cmd, logname, env or {}))
        return None

    def flush(self, stage):
        if not self.pending:
            return
        jobs, self.pending = self.pending, []
        sc = slurm_cfg(self.cfg)
        root = SLURM_DIR.resolve()
        (root / "logs").mkdir(parents=True, exist_ok=True)
        stamp = time.strftime("%Y%m%d-%H%M%S")
        size = int(sc["max_array_size"])
        chunks = [jobs[i:i + size] for i in range(0, len(jobs), size)]
        for c, chunk in enumerate(chunks):
            tag = f"{stage}-{stamp}" + (f"-{c + 1}" if len(chunks) > 1 else "")
            manifest = root / f"{tag}.jobs"
            scripts = [self._job_script(stage, manifest, *job) for job in chunk]
            manifest.write_text("".join(f"{s}\n" for s in scripts))
            script = root / f"{tag}.sbatch"
            script.write_text(self._array_script(stage, tag, manifest, len(chunk), sc))
            self._submit(script, manifest, chunk, sc)

    def _job_script(self, stage, manifest, wd, cmd, logname, env):
        name = Path(logname).stem
        q = shlex.quote
        exports = "".join(f'export {k}="${{{k}:-{v}}}"\n' for k, v in env.items())
        path = wd / f"{name}.sh"
        path.write_text(
            "#!/bin/bash\n"
            f"# {stage} job, task of {manifest.name}; written by pipeline.py\n"
            f"cd {q(str(wd))} || exit 1\n"
            f"rm -f {q(name + '.queued')}\n"
            f'echo "slurm ${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}} $(hostname)"'
            f" > {q(name + '.running')}\n"
            f"{exports}"
            f"{' '.join(q(str(a)) for a in cmd)} > {q(logname)} 2>&1\n"
            "rc=$?\n"
            f"rm -f {q(name + '.running')}\n"
            "exit $rc\n")
        path.chmod(0o755)
        return path

    def _array_script(self, stage, tag, manifest, n, sc):
        nprocs = self.cfg["nprocs"]
        mem = sc["mem_per_cpu_mb"] or int(1.25 * self.cfg["maxcore"])
        lines = ["#!/bin/bash",
                 f"#SBATCH --job-name=baobei-{stage}",
                 f"#SBATCH --array=1-{n}%{sc['max_parallel']}",
                 "#SBATCH --nodes=1"]
        if stage == "crest":     # one process, OpenMP threads
            lines += ["#SBATCH --ntasks=1", f"#SBATCH --cpus-per-task={nprocs}"]
        else:                    # ORCA %pal runs nprocs MPI ranks
            lines += [f"#SBATCH --ntasks={nprocs}", "#SBATCH --cpus-per-task=1"]
        lines += [f"#SBATCH --mem-per-cpu={mem}M",
                  f"#SBATCH --time={sc['time'].get(stage, '24:00:00')}",
                  f"#SBATCH --output={SLURM_DIR.resolve() / 'logs' / tag}_%a.log"]
        if sc["partition"]:
            lines.append(f"#SBATCH --partition={sc['partition']}")
        if sc["account"]:
            lines.append(f"#SBATCH --account={sc['account']}")
        lines += list(sc["sbatch_lines"])
        lines += ["", f"# {stage}: array task i runs the job script on line i of the manifest",
                  *sc["setup"]]
        if stage == "crest":
            lines.append(f"export OMP_NUM_THREADS={nprocs}")
        lines.append(f'exec bash "$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {shlex.quote(str(manifest))})"')
        return "\n".join(lines) + "\n"

    def _submit(self, script, manifest, chunk, sc):
        if not self.submit:
            print(f"[slurm] wrote {script} ({len(chunk)} jobs), not submitted")
            return
        markers = []
        for i, (wd, _, logname, _) in enumerate(chunk, 1):
            m = wd / f"{Path(logname).stem}.queued"
            m.write_text(f"queued {manifest} {i}\n")
            markers.append(m)
        try:
            p = subprocess.run([sc["sbatch"], "--parsable", str(script)],
                               capture_output=True, text=True)
        except OSError as exc:
            p = subprocess.CompletedProcess([], 1, "", str(exc))
        if p.returncode != 0:
            for m in markers:
                if m.exists():
                    m.unlink()
            print(f"[slurm] sbatch failed for {script}: {p.stderr.strip()}\n"
                  f"        submit it by hand with: sbatch {script}")
            return
        jobid = p.stdout.strip().split(";")[0]
        manifest.with_suffix(".jobid").write_text(jobid + "\n")
        print(f"[slurm] submitted {len(chunk)} jobs as array {jobid} ({script.name})")


RUNNER = Runner()


def orca_job(cfg, wd, name, keywords, xyz, v):
    """Write and run one ORCA job in wd unless it already finished."""
    wd = Path(wd)
    out = wd / f"{name}.out"
    if orca_done(out):
        return True
    if active(cfg, wd, name):
        print(f"  ORCA {name}: {wd} (already {orca_state(cfg, wd, name)})")
        return False
    inp = (f"! {keywords} CPCM({cfg['solvent_orca']}) TightSCF\n"
           f"%pal nprocs {cfg['nprocs']} end\n"
           f"%maxcore {cfg['maxcore']}\n"
           f"* xyzfile {v['charge']} {v['mult']} {xyz}\n")
    (wd / f"{name}.inp").write_text(inp)
    print(f"  ORCA {name}: {wd}")
    if RUNNER.add(wd, [cfg["orca"], f"{name}.inp"], f"{name}.out") is None:
        return False
    ok = orca_done(out)
    if not ok:
        print(f"  !! failed: {out}")
    return ok


# ---------------------------------------------------------------- selections
def sp_selection(cfg, xyzfile):
    """CREST conformers inside sp_window as [(index, energy, xyz_block)], and the total."""
    ens = read_multixyz(xyzfile)
    found = [e for e, _ in ens if e is not None]
    if len(found) < len(ens):
        print(f"  !! {xyzfile}: {len(ens) - len(found)} structures without an energy, ignored")
    if not found:
        return [], len(ens)
    emin = min(found)
    sel = [(k, e, block) for k, (e, block) in enumerate(ens)
           if e is not None and (e - emin) * H2KCAL <= cfg["sp_window"]]
    return sel, len(ens)


def opt_selection(cfg, sp, spec):
    """Finished single points pooled over variants, inside opt_window, capped at opt_max."""
    rows = []
    for var, v in spec["variants"].items():
        for d in sorted((WORK / sp / var / "sp").glob("c*")):
            if orca_done(d / "sp.out"):
                rows.append((final_energy(d / "sp.out"), var, d, v))
    if not rows:
        return [], 0
    rows.sort(key=lambda r: r[0])
    emin = rows[0][0]
    sel = [r for r in rows
           if (r[0] - emin) * H2KCAL <= cfg["opt_window"]][:cfg["opt_max"]]
    return sel, len(rows)


def best_conformer(spec, sp):
    """(G_Eh, variant, opt_dir) of the lowest screen-level G, or None."""
    best = None
    for var in spec["variants"]:
        for d in sorted((WORK / sp / var / "opt").glob("c*")):
            if not orca_done(d / "opt.out") or has_imaginary(d / "opt.out"):
                continue
            corr = g_minus_eel(d / "opt.out")
            if corr is None:
                continue
            g = final_energy(d / "opt.out") + corr
            if best is None or g < best[0]:
                best = (g, var, d)
    return best


# ---------------------------------------------------------------- stages
def stage_crest(cfg):
    for sp, spec, var, v in variants(cfg):
        wd = WORK / sp / var / "crest"
        wd.mkdir(parents=True, exist_ok=True)
        if crest_done(wd):
            print(f"[crest] skip {sp}/{var} (done)")
            continue
        if active(cfg, wd, "crest"):
            print(f"[crest] skip {sp}/{var} ({crest_state(cfg, wd)})")
            continue
        shutil.copy(v["xyz"], wd / "start.xyz")
        cmd = [cfg["crest"], "start.xyz", "--gfn2",
               "--alpb", cfg["solvent_xtb"],
               "--chrg", str(v["charge"]), "--uhf", str(v["mult"] - 1),
               "--T", str(cfg["nprocs"]), "--ewin", str(cfg["crest_ewin"])]
        if spec.get("nci"):
            cmd.append("--nci")
        print(f"[crest] {sp}/{var}")
        # conda-forge's pthreads OpenBLAS otherwise floods crest.log with warnings
        rc = RUNNER.add(wd, cmd, "crest.log", env={"OPENBLAS_NUM_THREADS": "1"})
        if rc not in (0, None):
            print(f"  !! CREST returned nonzero, check {wd/'crest.log'}")


def stage_sp(cfg):
    for sp, spec, var, v in variants(cfg):
        f = WORK / sp / var / "crest" / "crest_conformers.xyz"
        if active(cfg, f.parent, "crest"):
            print(f"[sp] CREST still running for {sp}/{var}, skipping")
            continue
        if not crest_done(f.parent):
            print(f"[sp] no finished CREST output for {sp}/{var}, skipping")
            continue
        sel, total = sp_selection(cfg, f)
        for k, e, block in sel:
            wd = WORK / sp / var / "sp" / f"c{k:03d}"
            wd.mkdir(parents=True, exist_ok=True)
            (wd / "conf.xyz").write_text(block)
            orca_job(cfg, wd, "sp", cfg["methods"]["screen"], "conf.xyz", v)
        print(f"[sp] {sp}/{var}: {len(sel)} of {total} conformers screened")


def stage_opt(cfg):
    for sp, spec in cfg["species"].items():
        busy = [d for var in spec["variants"]
                for d in (WORK / sp / var / "sp").glob("c*") if active(cfg, d, "sp")]
        if busy:
            print(f"[opt] {sp}: {len(busy)} single points still running or queued, skipping")
            continue
        sel, total = opt_selection(cfg, sp, spec)
        if not sel:
            print(f"[opt] {sp}: no finished single points")
            continue
        print(f"[opt] {sp}: optimizing {len(sel)} of {total} (pooled over variants)")
        for e, var, d, v in sel:
            wd = WORK / sp / var / "opt" / d.name
            wd.mkdir(parents=True, exist_ok=True)
            shutil.copy(d / "conf.xyz", wd / "conf.xyz")
            orca_job(cfg, wd, "opt", cfg["methods"]["screen"] + " Opt Freq",
                     "conf.xyz", v)


def stage_extras(cfg):
    """C4 descriptors on the lowest-G conformer of each complex and monomer."""
    for sp, spec in cfg["species"].items():
        busy = [d for var in spec["variants"]
                for d in (WORK / sp / var / "opt").glob("c*") if active(cfg, d, "opt")]
        if busy:
            print(f"[extras] {sp}: {len(busy)} optimizations still running or queued, skipping")
            continue
        best = best_conformer(spec, sp)
        if best is None:
            continue
        _, var, d = best
        v = spec["variants"][var]
        # oxidation potential: vertical radical cation at the neutral geometry (C4)
        if spec.get("oxidation"):
            wd = WORK / sp / var / "extras" / "cation"
            wd.mkdir(parents=True, exist_ok=True)
            shutil.copy(d / "opt.xyz", wd / "geom.xyz")
            vc = {"charge": v["charge"] + 1, "mult": v["mult"] + 1}
            orca_job(cfg, wd, "cation", cfg["methods"]["high"], "geom.xyz", vc)
        # counterpoise-corrected interaction energy (C4 energy decomposition input)
        if spec.get("fragments"):
            wd = WORK / sp / var / "extras" / "cp"
            wd.mkdir(parents=True, exist_ok=True)
            out = wd / "cp.out"
            if not orca_done(out) and not active(cfg, wd, "cp"):
                if not (wd / "geom.xyz").exists():   # keep hand-added fragment labels
                    shutil.copy(d / "opt.xyz", wd / "geom.xyz")
                inp = (f"! {cfg['methods']['high']} CPCM({cfg['solvent_orca']}) "
                       f"TightSCF LED\n"
                       f"%pal nprocs {cfg['nprocs']} end\n%maxcore {cfg['maxcore']}\n"
                       f"%method frag true end\n"
                       f"* xyzfile {v['charge']} {v['mult']} geom.xyz\n")
                (wd / "cp.inp").write_text(inp)
                print(f"  ORCA cp: {wd}  (needs fragment labels in geom.xyz)")
                RUNNER.add(wd, [cfg["orca"], "cp.inp"], "cp.out")


def stage_hl(cfg):
    for sp, spec, var, v in variants(cfg):
        for d in sorted((WORK / sp / var / "opt").glob("c*")):
            if not orca_done(d / "opt.out") or has_imaginary(d / "opt.out"):
                continue
            wd = WORK / sp / var / "hl" / d.name
            wd.mkdir(parents=True, exist_ok=True)
            shutil.copy(d / "opt.xyz", wd / "geom.xyz")
            orca_job(cfg, wd, "hl", cfg["methods"]["high"], "geom.xyz", v)


# ---------------------------------------------------------------- analysis
def species_ensemble(cfg, sp, spec, level):
    """Ensemble free energy G = -RT ln sum exp(-G_i/RT), plus per-conformer weights."""
    RT = R_KCAL * cfg["temperature"]
    confs, flagged, no_thermo = [], [], []
    for var in spec["variants"]:
        for d in sorted((WORK / sp / var / "opt").glob("c*")):
            out = d / "opt.out"
            if not orca_done(out):
                continue
            if has_imaginary(out):
                flagged.append(f"{var}/{d.name}")
                continue
            corr = g_minus_eel(out)
            if corr is None:
                no_thermo.append(f"{var}/{d.name}")
                continue
            if level == "high":
                hl = WORK / sp / var / "hl" / d.name / "hl.out"
                if not orca_done(hl):
                    continue
                e = final_energy(hl)
            else:
                e = final_energy(out)
            confs.append({"variant": var, "conf": d.name,
                          "G_kcal": (e + corr) * H2KCAL})
    if flagged:
        print(f"  {sp}: excluded (imaginary frequency): {', '.join(flagged)}")
    if no_thermo:
        print(f"  {sp}: excluded (no G-E(el) in opt.out): {', '.join(no_thermo)}")
    if not confs:
        return None, []
    # energy-based de-duplication of conformers that converged to the same minimum
    confs.sort(key=lambda c: c["G_kcal"])
    uniq = [confs[0]]
    for c in confs[1:]:
        if c["G_kcal"] - uniq[-1]["G_kcal"] > cfg["dedupe_kcal"]:
            uniq.append(c)
    gmin = uniq[0]["G_kcal"]
    z = [math.exp(-(c["G_kcal"] - gmin) / RT) for c in uniq]
    Z = sum(z)
    for c, zi in zip(uniq, z):
        c["dG_rel"] = c["G_kcal"] - gmin
        c["weight"] = zi / Z
    return gmin - RT * math.log(Z), uniq


def write_ranking(cfg, dGb):
    """One row per monomer: the computational side of C2 and C3.

    Fill measured_IF by hand once sensor data exists, then correlate.
    A monotonic fit supports C2; a peak in the middle supports C3.
    """
    rows = []
    for m in cfg.get("monomers", []):
        hist, hisd = f"{m} histamine", f"{m} histidine"
        if hist not in dGb:
            continue
        gap = dGb[hist] - dGb[hisd] if hisd in dGb else None
        ox = oxidation_potential(cfg, f"{m}_histamine_1to1")
        rows.append([m, f"{dGb[hist]:.2f}",
                     "" if gap is None else f"{gap:.2f}",
                     "" if ox is None else f"{ox:.2f}", ""])
    if not rows:
        return
    with open(RESULTS / "ranking.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["monomer", "dG_bind_histamine_kcal", "selectivity_gap_kcal",
                    "oxidation_potential_V_vs_SHE", "measured_IF"])
        w.writerows(rows)
    print(f"[analyze] wrote {RESULTS}/ranking.csv (C2/C3 table, measured_IF blank)")


def oxidation_potential(cfg, sp):
    """Vertical oxidation potential vs SHE, from the neutral and cation energies."""
    spec = cfg["species"].get(sp)
    if not spec:
        return None
    for var in spec["variants"]:
        cat = WORK / sp / var / "extras" / "cation" / "cation.out"
        neu = WORK / sp / var / "hl"
        if not orca_done(cat):
            continue
        geoms = [d / "hl.out" for d in sorted(neu.glob("c*")) if orca_done(d / "hl.out")]
        if not geoms:
            continue
        e_n = min(final_energy(g) for g in geoms)
        return (final_energy(cat) - e_n) * 27.2114 - cfg.get("she_ref_V", 4.28)
    return None


def stage_analyze(cfg):
    RESULTS.mkdir(exist_ok=True)
    T = cfg["temperature"]
    RT = R_KCAL * T
    ss = RT * math.log(0.082057366 * T)   # 1 atm -> 1 M, per species (~1.89 kcal/mol at 298 K)
    summary = []
    for level in ("screen", "high"):
        G = {}
        with open(RESULTS / f"conformers_{level}.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["species", "variant", "conf", "dG_rel_kcal", "weight"])
            for sp, spec in cfg["species"].items():
                g, confs = species_ensemble(cfg, sp, spec, level)
                if g is None:
                    continue
                G[sp] = g
                for c in confs:
                    w.writerow([sp, c["variant"], c["conf"],
                                f"{c['dG_rel']:.2f}", f"{c['weight']:.3f}"])
        dGb = {}
        for b in cfg["bindings"]:
            needed = [b["complex"], *b["parts"]]
            if not all(s in G for s in needed):
                continue
            n_parts = sum(b["parts"].values())
            dg = G[b["complex"]] - sum(n * G[s] for s, n in b["parts"].items())
            dg += ss * (1 - n_parts)
            dGb[b["name"]] = dg
            summary.append([level, "binding", b["name"], f"{dg:.2f}", ""])
        for a, bname in cfg["selectivity"]:
            if a in dGb and bname in dGb:
                ddg = dGb[a] - dGb[bname]
                summary.append([level, "selectivity gap (C1)", f"{a} vs {bname}",
                                f"{ddg:.2f}", f"{math.exp(-ddg / RT):.2f}"])
        if level == "high":
            write_ranking(cfg, dGb)
    with open(RESULTS / "summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["level", "type", "name", "dG_kcal_1M", "K_ratio"])
        w.writerows(summary)
    for row in summary:
        print("  " + " | ".join(row))
    print(f"[analyze] wrote {RESULTS}/summary.csv and conformers_*.csv")


# ---------------------------------------------------------------- status
STATES = ("done", "running", "queued", "failed", "not started")


def _job_dirs(parent, name):
    """Subdirectories of parent that hold any trace of job `name`."""
    exts = (".inp", ".out", ".running", ".queued")
    return {d.name for d in Path(parent).glob("c*")
            if any((d / f"{name}{x}").exists() for x in exts)}


def species_status(cfg, sp, spec):
    """{stage: Counter(state)} and the paths of failed jobs, without touching anything."""
    counts = {s: Counter() for s in ("crest", "sp", "opt", "hl", "extras")}
    failed = []

    def add(stage, state, path):
        counts[stage][state] += 1
        if state == "failed":
            failed.append(path)

    for var in spec["variants"]:
        base = WORK / sp / var
        cw = base / "crest"
        add("crest", crest_state(cfg, cw), cw / "crest.log")
        f = cw / "crest_conformers.xyz"
        names = _job_dirs(base / "sp", "sp")
        if crest_done(cw):
            names |= {f"c{k:03d}" for k, _, _ in sp_selection(cfg, f)[0]}
        for n in sorted(names):
            add("sp", orca_state(cfg, base / "sp" / n, "sp"), base / "sp" / n / "sp.out")
    want = {(var, d.name) for _, var, d, _ in opt_selection(cfg, sp, spec)[0]}
    for var in spec["variants"]:
        want |= {(var, n) for n in _job_dirs(WORK / sp / var / "opt", "opt")}
    for var, n in sorted(want):
        d = WORK / sp / var / "opt" / n
        add("opt", orca_state(cfg, d, "opt"), d / "opt.out")
        if orca_done(d / "opt.out") and not has_imaginary(d / "opt.out"):
            h = WORK / sp / var / "hl" / n
            add("hl", orca_state(cfg, h, "hl"), h / "hl.out")
    extras = set()
    best = best_conformer(spec, sp)
    if best is not None:
        if spec.get("oxidation"):
            extras.add((best[1], "cation"))
        if spec.get("fragments"):
            extras.add((best[1], "cp"))
    for var in spec["variants"]:
        for name in ("cation", "cp"):
            wd = WORK / sp / var / "extras" / name
            if any((wd / f"{name}{x}").exists() for x in (".inp", ".out", ".running", ".queued")):
                extras.add((var, name))
    for var, name in sorted(extras):
        wd = WORK / sp / var / "extras" / name
        add("extras", orca_state(cfg, wd, name), wd / f"{name}.out")
    return counts, failed


def stage_status(cfg):
    width = max(len("species"), *(len(s) for s in cfg["species"])) + 2
    cols = [max(len(s), 4) + 2 for s in STATES]
    head = f"{'species':<{width}}{'stage':<8}" + "".join(f"{s:>{n}}" for s, n in zip(STATES, cols))
    print(head)
    print("-" * len(head))
    total = {s: Counter() for s in ("crest", "sp", "opt", "hl", "extras")}
    all_failed = []
    for sp, spec in cfg["species"].items():
        counts, failed = species_status(cfg, sp, spec)
        all_failed += failed
        label = sp
        for stage, c in counts.items():
            total[stage].update(c)
            if not sum(c.values()):
                continue
            print(f"{label:<{width}}{stage:<8}" + "".join(f"{c[s]:>{n}}" for s, n in zip(STATES, cols)))
            label = ""
    print("-" * len(head))
    for stage, c in total.items():
        print(f"{'TOTAL':<{width}}{stage:<8}" + "".join(f"{c[s]:>{n}}" for s, n in zip(STATES, cols)))
    if all_failed:
        print(f"\nfailed jobs ({len(all_failed)}):")
        for p in all_failed[:30]:
            print(f"  {p}")
        if len(all_failed) > 30:
            print(f"  ... and {len(all_failed) - 30} more")
    if squeue_ids(cfg) is None and any(c["running"] + c["queued"] for c in total.values()):
        print("\nnote: squeue not available, so SLURM jobs that died without cleaning up "
              "still count as running or queued")


STAGES = {"crest": stage_crest, "sp": stage_sp, "opt": stage_opt,
          "hl": stage_hl, "extras": stage_extras, "analyze": stage_analyze}


def main(argv=None):
    global RUNNER
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=[*STAGES, "all", "status"])
    ap.add_argument("config", nargs="?", default="config.json")
    ap.add_argument("--slurm", action="store_true",
                    help="write SLURM job arrays and submit them instead of running jobs here")
    ap.add_argument("--no-submit", action="store_true",
                    help="with --slurm: write the scripts under slurm/ but do not call sbatch")
    args = ap.parse_args(argv)
    if args.stage == "all" and args.slurm:
        ap.error("--slurm runs one stage at a time: each stage needs the previous one "
                 "finished. Submit a stage, wait for `status`, then submit the next.")
    cfg = load_cfg(args.config)
    if args.stage == "status":
        stage_status(cfg)
        return
    RUNNER = Runner(cfg, slurm=args.slurm, submit=not args.no_submit)
    todo = list(STAGES) if args.stage == "all" else [args.stage]
    for s in todo:
        try:
            STAGES[s](cfg)
        except FileNotFoundError as exc:
            if exc.filename and Path(str(exc.filename)).name in (
                    Path(cfg["orca"]).name, Path(cfg["crest"]).name):
                sys.exit(f"cannot run {exc.filename}: set the full path in {args.config}")
            raise
        if args.slurm:
            RUNNER.flush(s)


if __name__ == "__main__":
    main()
