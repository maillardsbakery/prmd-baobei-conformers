#!/usr/bin/env python3
"""Write placeholder start geometries for every variant in a config.

  python tests/make_fake_inputs.py config_water.json <run_dir>

Only for dry runs against tests/fake_bin: the atoms have the right elements
and counts (so the fake programs give stoichiometrically consistent energies),
but the coordinates are a grid, not a molecule. MAA gets the real
GFN2-xTB geometry from tests/fixtures/crest3/maa_water_full/start.xyz.
Real inputs come from Avogadro, as the README says.
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MONOMERS = {"MAA": "C4H6O2", "AA": "C3H4O2", "IA": "C5H6O4", "AM": "C3H5NO",
            "4VP": "C7H7N", "2VP": "C7H7N", "AMPS": "C7H13NO4S", "HEMA": "C6H10O3"}
TEMPLATES = {"histamine": {"neutral_tau": "C5H9N3", "neutral_pi": "C5H9N3",
                           "cation": "C5H10N3"},
             "histidine": {"zwitterion": "C6H9N3O2", "cation": "C6H10N3O2"}}
NEUTRAL = {"histamine": "neutral_tau", "histidine": "zwitterion"}


def atoms(formula):
    out = []
    for el, k in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        out += [el] * int(k or 1)
    return out


def grid_xyz(elements, title):
    rows = [f"{el:<2} {1.5 * (i % 4):10.4f} {1.5 * (i // 4 % 4):10.4f} {1.5 * (i // 16):10.4f}"
            for i, el in enumerate(elements)]
    return f"{len(elements)}\n{title}\n" + "\n".join(rows) + "\n"


def formula_for(species, variant):
    if species in MONOMERS:
        return atoms(MONOMERS[species])
    if species in TEMPLATES:
        return atoms(TEMPLATES[species][variant])
    m = re.fullmatch(r"(.+)_(histamine|histidine)_1to1", species)
    if not m:
        raise SystemExit(f"don't know the composition of {species}")
    mono, tpl = m.groups()
    tv = variant[len("hbond_"):] if variant.startswith("hbond_") else NEUTRAL[tpl]
    return atoms(MONOMERS.get(mono, "C4H6O2")) + atoms(TEMPLATES[tpl][tv])


def main(cfg_path, run_dir):
    cfg = json.loads(Path(cfg_path).read_text())
    run_dir = Path(run_dir)
    for sp, spec in cfg["species"].items():
        for var, v in spec["variants"].items():
            path = run_dir / v["xyz"]
            path.parent.mkdir(parents=True, exist_ok=True)
            if sp == "MAA":
                path.write_text((HERE / "fixtures/crest3/maa_water_full/start.xyz").read_text())
            else:
                path.write_text(grid_xyz(formula_for(sp, var), f"FAKE {sp}/{var}"))


if __name__ == "__main__":
    main(*sys.argv[1:3])
