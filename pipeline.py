#!/usr/bin/env python3
"""
Baobei 2.0 conformer pipeline

  crest    : CREST (GFN2-xTB, ALPB solvent) conformer search per species/variant
  sp       : r2SCAN-3c single points on every CREST conformer inside sp_window
  opt      : pool all variants of a species, rank by SP energy, Opt+Freq the best
  hl       : high-level single points (second functional) on optimized geometries
  analyze  : ensemble free energies, binding free energies, selectivity
  all      : run every stage in order

Every stage is resumable: finished jobs (ORCA TERMINATED NORMALLY / CREST output
present) are skipped, so you can stop, add compute, and rerun.

Usage:  python pipeline.py <stage> [config.json]
"""
import csv
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

H2KCAL = 627.509474          # Eh -> kcal/mol
R_KCAL = 0.0019872043        # kcal/(mol K)
WORK = Path("work")
RESULTS = Path("results")


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


def run(cmd, cwd, logname):
    with open(Path(cwd) / logname, "w") as log:
        return subprocess.run(cmd, cwd=cwd, stdout=log,
                              stderr=subprocess.STDOUT).returncode


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


def orca_job(cfg, wd, name, keywords, xyz, v):
    """Write and run one ORCA job in wd unless it already finished."""
    wd = Path(wd)
    out = wd / f"{name}.out"
    if orca_done(out):
        return True
    inp = (f"! {keywords} CPCM({cfg['solvent_orca']}) TightSCF\n"
           f"%pal nprocs {cfg['nprocs']} end\n"
           f"%maxcore {cfg['maxcore']}\n"
           f"* xyzfile {v['charge']} {v['mult']} {xyz}\n")
    (wd / f"{name}.inp").write_text(inp)
    print(f"  ORCA {name}: {wd}")
    run([cfg["orca"], f"{name}.inp"], wd, f"{name}.out")
    ok = orca_done(out)
    if not ok:
        print(f"  !! failed: {out}")
    return ok


# ---------------------------------------------------------------- stages
def stage_crest(cfg):
    for sp, spec, var, v in variants(cfg):
        wd = WORK / sp / var / "crest"
        wd.mkdir(parents=True, exist_ok=True)
        if (wd / "crest_conformers.xyz").exists():
            print(f"[crest] skip {sp}/{var} (done)")
            continue
        shutil.copy(v["xyz"], wd / "start.xyz")
        cmd = [cfg["crest"], "start.xyz", "--gfn2",
               "--alpb", cfg["solvent_xtb"],
               "--chrg", str(v["charge"]), "--uhf", str(v["mult"] - 1),
               "--T", str(cfg["nprocs"]), "--ewin", str(cfg["crest_ewin"])]
        if spec.get("nci"):
            cmd.append("--nci")
        print(f"[crest] {sp}/{var}")
        if run(cmd, wd, "crest.log") != 0:
            print(f"  !! CREST returned nonzero, check {wd/'crest.log'}")


def stage_sp(cfg):
    for sp, spec, var, v in variants(cfg):
        f = WORK / sp / var / "crest" / "crest_conformers.xyz"
        if not f.exists():
            print(f"[sp] no CREST output for {sp}/{var}, skipping")
            continue
        ens = read_multixyz(f)
        emin = min(e for e, _ in ens)
        kept = 0
        for k, (e, block) in enumerate(ens):
            if (e - emin) * H2KCAL > cfg["sp_window"]:
                continue
            wd = WORK / sp / var / "sp" / f"c{k:03d}"
            wd.mkdir(parents=True, exist_ok=True)
            (wd / "conf.xyz").write_text(block)
            orca_job(cfg, wd, "sp", cfg["methods"]["screen"], "conf.xyz", v)
            kept += 1
        print(f"[sp] {sp}/{var}: {kept} of {len(ens)} conformers screened")


def stage_opt(cfg):
    for sp, spec in cfg["species"].items():
        rows = []
        for var, v in spec["variants"].items():
            for d in sorted((WORK / sp / var / "sp").glob("c*")):
                if orca_done(d / "sp.out"):
                    rows.append((final_energy(d / "sp.out"), var, d, v))
        if not rows:
            print(f"[opt] {sp}: no finished single points")
            continue
        rows.sort(key=lambda r: r[0])
        emin = rows[0][0]
        sel = [r for r in rows
               if (r[0] - emin) * H2KCAL <= cfg["opt_window"]][:cfg["opt_max"]]
        print(f"[opt] {sp}: optimizing {len(sel)} of {len(rows)} (pooled over variants)")
        for e, var, d, v in sel:
            wd = WORK / sp / var / "opt" / d.name
            wd.mkdir(parents=True, exist_ok=True)
            shutil.copy(d / "conf.xyz", wd / "conf.xyz")
            orca_job(cfg, wd, "opt", cfg["methods"]["screen"] + " Opt Freq",
                     "conf.xyz", v)


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
    confs, flagged = [], []
    for var in spec["variants"]:
        for d in sorted((WORK / sp / var / "opt").glob("c*")):
            out = d / "opt.out"
            if not orca_done(out):
                continue
            if has_imaginary(out):
                flagged.append(f"{var}/{d.name}")
                continue
            corr = g_minus_eel(out)
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
                summary.append([level, "selectivity", f"{a} vs {bname}",
                                f"{ddg:.2f}", f"{math.exp(-ddg / RT):.2f}"])
    with open(RESULTS / "summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["level", "type", "name", "dG_kcal_1M", "K_ratio"])
        w.writerows(summary)
    for row in summary:
        print("  " + " | ".join(row))
    print(f"[analyze] wrote {RESULTS}/summary.csv and conformers_*.csv")


STAGES = {"crest": stage_crest, "sp": stage_sp, "opt": stage_opt,
          "hl": stage_hl, "analyze": stage_analyze}

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in [*STAGES, "all"]:
        sys.exit(__doc__)
    cfg = load_cfg(sys.argv[2] if len(sys.argv) > 2 else "config.json")
    todo = list(STAGES) if sys.argv[1] == "all" else [sys.argv[1]]
    for s in todo:
        STAGES[s](cfg)
