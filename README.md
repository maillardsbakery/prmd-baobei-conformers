# Baobei 2.0 conformer pipeline

CREST (GFN2-xTB) explores, r2SCAN-3c screens and optimizes, wB97M-V re-scores, and the
analysis turns everything into ensemble binding free energies and a histamine vs
histidine selectivity.

## What it does, per species

1. **crest**: conformer search for each variant (tautomer, protonation state) in implicit
   porogen. Complexes run with `--nci` so the partners stay together.
2. **sp**: r2SCAN-3c single point on every CREST conformer within `sp_window` (10 kcal/mol).
   This is the "don't trust xTB's ranking" safety net.
3. **opt**: all variants of a species are pooled and ranked together, so a neutral pair and
   an ion pair compete on equal footing. The best `opt_max` within `opt_window` get
   r2SCAN-3c Opt + Freq.
4. **hl**: wB97M-V/def2-TZVP single point on each optimized geometry.
5. **analyze**: for both levels (screen and high), computes
   - G_i = E_el + (G - E_el) from the r2SCAN-3c frequencies
   - ensemble free energy G = -RT ln sum exp(-G_i/RT)
   - binding: dG = G(complex) - sum n G(part), corrected from 1 atm to 1 M
   - selectivity: ddG = dG(histamine) - dG(histidine), and K ratio = exp(-ddG/RT)

Reporting both levels side by side is the functional sensitivity check: if the
selectivity sign flips between r2SCAN-3c and wB97M-V, the result is not robust yet.

## Setup

- CREST 3 and xtb (conda: `conda install -c conda-forge crest xtb`)
- ORCA 6 (free academic licence). Set `"orca"` to the **full path**, which ORCA needs for
  parallel runs, and install OpenMPI to match.
- Python 3.8+, no extra packages.

Build starting geometries in Avogadro and save them under `inputs/<species>/<variant>.xyz`
to match `config.json`. They only need to be sensible, not optimized. For ion pairs, move
the carboxylic proton onto the histamine amine before saving.

## Run

```
python pipeline.py crest
python pipeline.py sp
python pipeline.py opt
python pipeline.py hl
python pipeline.py analyze
```

or `python pipeline.py all`. Every stage skips finished jobs, so you can interrupt and
rerun. Set `solvent_xtb` and `solvent_orca` to your porogen once it is chosen.

## Outputs

- `results/summary.csv`: binding and selectivity at both levels
- `results/conformers_screen.csv`, `conformers_high.csv`: every conformer kept, its
  relative free energy and Boltzmann weight, and which variant it came from

## Things the script will not catch for you

- **MAA self-association.** The reference state here is n separate MAA molecules. MAA
  dimerizes strongly in low-polarity porogens, which would make real binding weaker than
  reported. To test it, add an `MAA_dimer` species and use it as a part in `bindings`.
- **Ion pairs in implicit solvent.** Their stability is sensitive to the solvation model.
  If ion pairs dominate, confirm with one or two explicit porogen molecules.
- **De-duplication is energy-only** (0.05 kcal/mol). Check the top structures by eye.
- **Imaginary frequencies.** Conformers with one are excluded and listed in the
  analyze output. Re-optimize them from a displaced geometry if they were low in energy.
- Low-frequency modes inflate entropy. Check how your ORCA version treats them
  (quasi-RRHO) and state it in the methods.

## Rough cost

Histamine + 4 MAA is about 65 atoms. Expect the CREST searches to take hours on 8 cores,
and the DFT stages to dominate: a few hundred r2SCAN-3c single points plus up to 30
Opt + Freq jobs per species. Variants multiply this linearly, which is where bought
compute pays off.
