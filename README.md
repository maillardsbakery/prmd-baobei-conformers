# prmd-baobei-conformers

The computation column of the thesis: CREST (GFN2-xTB) explores, r2SCAN-3c screens and
optimizes, wB97M-V re-scores, and the analysis produces the numbers C1 to C4 are about.

| Hypothesis | What this repo produces | What it still needs |
| --- | --- | --- |
| **C1** selectivity gap above zero in implicit water | `results/summary.csv`, one gap per monomer | nothing; this repo can settle it |
| **C2** binding free energy predicts sensor performance in rank order | `results/ranking.csv`, monomers ranked by dG | measured imprinting factors, pasted into `measured_IF` |
| **C3** or performance peaks at a middle binding strength | same table; C2 and C3 are read off the same correlation | same |
| **C4** other descriptors beat binding energy alone | `selectivity_gap_kcal`, `oxidation_potential_V_vs_SHE`, LED output under `work/*/extras/cp/` | same |

C2 and C3 cannot be settled here. This repo produces the x-axis; the sensor work
produces the y-axis. Keep `measured_IF` blank until the electrodes exist.

## Solvent: run both

C1 is stated in implicit water, which is where the sensor works, so `config_water.json`
is the primary run. But the polymer forms in porogen, so `config_porogen.json` runs the
same species in acetonitrile. If the two disagree on monomer ranking, report both and say
which solvent each claim refers to. That disagreement is a finding, not a bug.

Protonation matters more in water than in porogen. At pH 7 histamine is largely
protonated, so the cation variant is included for both templates and competes with the
neutral forms in the same ranking.

## Run

```
python make_config.py water   > config_water.json
python make_config.py porogen > config_porogen.json

python pipeline.py crest   config_water.json
python pipeline.py sp      config_water.json
python pipeline.py opt     config_water.json
python pipeline.py hl      config_water.json
python pipeline.py extras  config_water.json
python pipeline.py analyze config_water.json
```

Every stage skips finished jobs, so you can interrupt and rerun. Edit `MONOMERS` in
`make_config.py` to change the set.

## Cost warning

The water config is 26 species: 8 monomers, 2 templates with their variants, and 16
one-to-one complexes. That is roughly 60 CREST searches and several hundred DFT jobs per
solvent, and two solvents doubles it. This is cluster work, not laptop work.

Start with three monomers spanning weak to strong binding (for example HEMA, MAA, AMPS)
and confirm the pipeline end to end before expanding. Three points cannot distinguish C2
from C3; you need at least five or six monomers for the shape of the curve to mean
anything.

## Setup

- CREST 3 and xtb (`conda install -c conda-forge crest xtb`)
- ORCA 6 (free academic licence). Set `"orca"` to the full path.
- Python 3.9+ (the `|` dict merge in `make_config.py` needs 3.9)

Build starting geometries in Avogadro under `inputs/<species>/<variant>.xyz`. For ion
pairs, move the acidic proton onto the histamine amine first. The LED fragment job needs
fragment labels in the xyz; the script writes the input but you must add those labels.

## Things the script will not catch

- **Monomer self-association.** The reference state is separate monomers. Carboxylic acids
  dimerize in porogen, which would make real binding weaker than reported. Add a dimer
  species if the porogen run matters to a claim.
- **Oxidation potentials are vertical**, computed at the neutral geometry with no
  relaxation, and the SHE reference (4.28 V) is a choice. Treat the numbers as a ranking,
  not as absolute potentials.
- **De-duplication is energy-only** (0.05 kcal/mol). Check top structures by eye.
- **Imaginary frequencies.** Those conformers are excluded and listed by `analyze`.
- Low-frequency modes inflate entropy; check how your ORCA version handles quasi-RRHO and
  state it in the methods.

## Untested

The parsers have never seen real CREST or ORCA output. The first cloud session's job is to
build fixtures and test them. Expect to fix a regex.
