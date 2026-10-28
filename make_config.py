#!/usr/bin/env python3
"""
Write a config covering every monomer, so C2 and C3 can be tested.

  python make_config.py water   > config_water.json    # C1, C2, C3, C4
  python make_config.py porogen > config_porogen.json  # prepolymerization check

C1 is stated in implicit water, which is where the sensor works. The polymer is
actually formed in porogen, so run the porogen config too and report both.
If the two disagree on monomer ranking, that disagreement is a result, not an
error: say which solvent each claim refers to.
"""
import json
import sys

# Monomers to rank. Keep the keys as directory-safe slugs.
MONOMERS = {
    "MAA":   {"name": "methacrylic acid",            "charge": 0},
    "AA":    {"name": "acrylic acid",                "charge": 0},
    "IA":    {"name": "itaconic acid",               "charge": 0},
    "AM":    {"name": "acrylamide",                  "charge": 0},
    "4VP":   {"name": "4-vinylpyridine",             "charge": 0},
    "2VP":   {"name": "2-vinylpyridine",             "charge": 0},
    "AMPS":  {"name": "2-acrylamido-2-methylpropanesulfonic acid", "charge": 0},
    "HEMA":  {"name": "2-hydroxyethyl methacrylate", "charge": 0},
}

# Templates. In water at pH 7 histamine is predominantly protonated, so the
# cation is the primary variant and the neutral is kept as a check.
TEMPLATES = {
    "histamine": {
        "neutral_tau":  {"charge": 0, "mult": 1},
        "neutral_pi":   {"charge": 0, "mult": 1},
        "cation":       {"charge": 1, "mult": 1},
    },
    "histidine": {
        "zwitterion":   {"charge": 0, "mult": 1},
        "cation":       {"charge": 1, "mult": 1},
    },
}

SOLVENTS = {
    "water":   {"xtb": "water",        "orca": "water"},
    "porogen": {"xtb": "acetonitrile", "orca": "acetonitrile"},
}


def build(solvent):
    species, bindings, selectivity = {}, [], []

    for m, meta in MONOMERS.items():
        species[m] = {
            "nci": False,
            "note": meta["name"],
            "oxidation": True,          # C4: oxidation potential of the monomer
            "variants": {
                "neutral": {"xyz": f"inputs/{m}/neutral.xyz",
                            "charge": meta["charge"], "mult": 1}
            },
        }

    for t, variants in TEMPLATES.items():
        species[t] = {
            "nci": False,
            "oxidation": t == "histamine",
            "variants": {v: {"xyz": f"inputs/{t}/{v}.xyz", **c}
                         for v, c in variants.items()},
        }

    for m in MONOMERS:
        for t, variants in TEMPLATES.items():
            key = f"{m}_{t}_1to1"
            species[key] = {
                "nci": True,
                "fragments": True,       # C4: counterpoise / LED decomposition
                "oxidation": t == "histamine",
                "variants": {
                    f"hbond_{v}": {"xyz": f"inputs/{key}/hbond_{v}.xyz", **c}
                    for v, c in variants.items()
                } | {
                    "ionpair": {"xyz": f"inputs/{key}/ionpair.xyz",
                                "charge": 0, "mult": 1}
                },
            }
            bindings.append({"name": f"{m} {t}", "complex": key,
                             "parts": {t: 1, m: 1}})
        selectivity.append([f"{m} histamine", f"{m} histidine"])

    return {
        "crest": "crest",
        "orca": "/opt/orca/orca",
        "nprocs": 8,
        "maxcore": 3000,
        "temperature": 298.15,
        "solvent_xtb": SOLVENTS[solvent]["xtb"],
        "solvent_orca": SOLVENTS[solvent]["orca"],
        "she_ref_V": 4.28,
        "crest_ewin": 10.0,
        "sp_window": 10.0,
        "opt_window": 4.0,
        "opt_max": 30,
        "dedupe_kcal": 0.05,
        "methods": {
            "screen": "r2SCAN-3c",
            "high": "wB97M-V def2-TZVP def2/J RIJCOSX",
        },
        "monomers": list(MONOMERS),
        "species": species,
        "bindings": bindings,
        "selectivity": selectivity,
    }


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "water"
    if which not in SOLVENTS:
        sys.exit(f"pick one of: {', '.join(SOLVENTS)}")
    print(json.dumps(build(which), indent=2))
