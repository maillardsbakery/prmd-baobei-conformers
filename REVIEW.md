# Review: first test pass of pipeline.py

What this session did, what it could not check, and what looks wrong. Chemistry
choices (functionals, solvent, windows, monomer set, protonation states) are
**unchanged**. Where I think one is wrong, it is in "Things I think are wrong"
below, not in the code.

## Summary

| Task | Result |
| --- | --- |
| 1. Parser tests | 44 tests in `tests/test_parsers.py`, all passing (56 in the whole suite). CREST fixtures are **genuine** CREST 3.0.2 output. ORCA fixtures are **hand-made** because ORCA needs a licensed download. |
| 2. Fake executables, end-to-end run | `tests/fake_bin/{crest,orca}` emit the fixtures. `python pipeline.py all config_water.json` runs end to end (about 25 s). `summary.csv` has 48 rows (2 levels × 16 bindings + 8 gaps) and `ranking.csv` has 8 monomers × 5 columns with `measured_IF` blank. Those shapes are what the config implies. |
| 3. `--slurm` | Writes one job array per stage under `slurm/` and submits it with `sbatch`. Tested against fake `sbatch`/`squeue`, and once with real CREST running inside the generated job script. **Not tested on a real SLURM cluster.** |
| 4. `status` | `python pipeline.py status config.json` prints, per species and stage, how many jobs are done, running, queued, failed or not started, then lists the failed jobs. |
| 5. Real CREST | `conda install -c conda-forge crest xtb` worked (via Miniforge; the container had no conda). CREST 3.0.2 / xtb 6.7.1 ran on MAA in water, both `--quick` and a full run. `read_multixyz` reads the genuine output correctly. Two real problems turned up; see "Found with real CREST". |

Run the tests with `pytest tests/` (about 70 s). The suite also passes on Python 3.9.

## Things I think are wrong (not changed)

### 1. Variants with different formulas are pooled into one ensemble. This matters most.

`stage_opt` ranks all variants of a species by raw single-point energy and keeps
those within `opt_window` (4 kcal/mol) of the lowest. `species_ensemble` then
Boltzmann-weights all of them together. That only makes sense when all variants
have the same atoms and charge. They do not:

| Species | Variants | Formulas |
| --- | --- | --- |
| histamine | neutral_tau, neutral_pi / cation | C5H9N3 (0) / C5H10N3 (+1) |
| histidine | zwitterion / cation | C6H9N3O2 (0) / C6H10N3O2 (+1) |
| X_histamine_1to1 | hbond_neutral_tau, hbond_neutral_pi, ionpair / hbond_cation | X+C5H9N3 (0) / X+C5H10N3 (+1) |
| X_histidine_1to1 | hbond_zwitterion, ionpair / hbond_cation | X+C6H9N3O2 (0) / X+C6H10N3O2 (+1) |

The protonated form has one more proton, so its total electronic energy is lower
by roughly the proton affinity plus solvation: hundreds of kcal/mol, not a few.
As a result:

- No neutral, zwitterion or ion-pair conformer is ever optimized. All
  `opt_max` slots go to the +1 variant. The fake run shows this: `work/histamine/neutral_*/opt/`
  and `work/*_1to1/ionpair/opt/` stay empty, and every mixed species has weight
  1.000 on its cation variant in `conformers_*.csv`.
- The README says the cation "competes with the neutral forms in the same ranking".
  It does not compete. It always wins, for an arithmetic reason, not a chemical one.
- C1 ends up comparing MAA·histamineH⁺ against MAA·histidineH⁺ (the histidine
  **cation**). At pH 7, histidine is mostly the zwitterion; with an imidazole pKa
  of about 6, only about a tenth is protonated. Histamine as the monocation is reasonable at pH 7.
- The binding stoichiometry still balances, but only because both sides happen to
  pick the +1 form. If one side ever lost its cation variant (for example a failed
  CREST run), ΔG_bind would silently be off by a proton.

Possible fixes (the choice is yours):

- (a) Make each protonation state its own species and combine them through
  experimental pKa values.
- (b) Keep them as variants, but only pool variants with the same formula and
  charge. Then pick or weight the states by pKa at pH 7.
- (c) Use an explicit reference free energy for the proton in water.

Whichever you choose, a check that every pooled variant has the same formula
would catch this class of error. The fake inputs (`tests/make_fake_inputs.py`)
already encode the formulas, if you want to start from there.

### 2. Protonation state of the monomers in water

Only neutral monomers are included. In water at pH 7:

- MAA (pKa ≈ 4.7) and AA are mostly carboxylate.
- AMPS is a sulfonic acid and fully deprotonated.
- IA has two acid groups.

The `ionpair` complex variants exist, but their reference state is the neutral
monomer plus histamine in the same ensemble. Item 1 means they never contribute
anyway. This may be deliberate for prepolymerization in porogen. For the water
claim (C1), it deserves a sentence in the methods.

### 3. Sign of the selectivity gap vs. the wording of C1

The code computes `gap = ΔG_bind(histamine) − ΔG_bind(histidine)`, with
`K_ratio = exp(−gap/RT) = K_histamine / K_histidine`. Under that definition:

- A **negative** gap means the monomer prefers histamine.
- C1, as written in the README and the docstring, is "selectivity gap **above zero**".

Either C1 means ΔG(histidine) − ΔG(histamine) > 0, and the code's sign is
flipped, or the claim should read "below zero". Please confirm before reading
C1 off `summary.csv`.

### 4. The oxidation potential is not strictly vertical

`stage_extras` puts the radical cation on the **lowest screen-level G**
conformer (r2SCAN-3c G). `oxidation_potential` then subtracts the **lowest wB97M-V
energy among all `hl` conformers** of that variant. That can be a different geometry,
so E(cation, geometry A) − E(neutral, geometry B) is not a vertical IP. A minimal
fix: use `hl/<same conformer>/hl.out` for the neutral. Also, if a species' best
variant is already a cation (histamine, and every histamine complex given item 1),
the "oxidation" makes a +2 doublet.

### 5. `ranking.csv` reports the complex's oxidation potential, not the monomer's

`write_ranking` calls `oxidation_potential(cfg, f"{m}_histamine_1to1")`. All 8
monomers have `"oxidation": true` and get a cation job, but their potentials
are never written anywhere. Which one does C4 need? The monomer's is the usual
descriptor for electropolymerization and stability. The complex's is what this
column holds now.

### 6. The "counterpoise" job is not counterpoise, and probably not valid input

- `stage_extras` labels the `cp` job "counterpoise-corrected interaction
  energy", but the input asks for `LED`, with no ghost-atom calculations.
- I am not confident ORCA 6 accepts `LED` with a DFT functional (it is documented for
  DLPNO-CCSD(T) and HFLD; there may be a DFT/HF variant I don't know).
- I don't recognise `%method frag true end` as ORCA input; fragments are usually
  assigned with `%geom Frags ... end end` or `(n)` labels on the atoms.
- Nothing parses `cp.out`, so C4 currently gets no number from this job.

Run one by hand before relying on it.

### 7. Smaller points

- `sp_window` (10) equals `crest_ewin` (10). CREST has already applied that window
  to the same GFN2 energies, so `sp_window` never removes anything on real output.
  Not wrong, but it is not doing the screening the docstring suggests.
- De-duplication compares each conformer only to the last kept one. Combined with
  the energy-only criterion the README already flags, two distinct conformers 0.04 kcal/mol
  apart collapse into one.
- The 4.28 V SHE reference and the 27.2114 eV/Eh conversion are fine. The 1 atm → 1 M
  correction (`ss * (1 − n_parts)`, ≈ −1.89 kcal/mol for 1:1) has the right sign.

## Found with real CREST (task 5)

Setup: Miniforge was installed and `conda create -c conda-forge crest xtb`
resolved crest 3.0.2, xtb 6.7.1 and tblite 0.6.0. MAA was built by hand,
pre-optimized with `xtb --opt --gfn 2 --alpb water`, and run through CREST:

| Run | Time | Conformers kept | `read_multixyz` |
| --- | --- | --- | --- |
| `--quick` before `--ewin 10` | 34 s on 2 threads | 4 | correct |
| `--quick` after `--ewin 10` (window became 5) | 52 s on 4 threads | 4 | correct |
| full, via `pipeline.py crest` | 76 s on 4 threads | 2 | correct |

`read_multixyz` also reads `crest_best.xyz` and `crest_rotamers.xyz`, and its
relative energies match CREST's own `crest.energies` to 0.002 kcal/mol.

1. **CREST 3 writes `crest_conformers.xyz` before it finishes.** It rewrites the
   file after every cycle. A snapshot taken mid-run had 3 structures; the finished
   run had 2. The pipeline treated "file exists" as done. So a CREST run killed by
   wall time would have been skipped forever, and `sp` would have screened a partial
   ensemble. **Fixed:** a CREST job now counts as done only if `crest.log` also
   contains `CREST terminated normally.` The test uses the genuine mid-run snapshot
   (`tests/fixtures/crest3/maa_water_midrun/`).
2. **The OpenBLAS warning flood.** conda-forge's pthreads OpenBLAS inside OpenMP CREST
   printed 184,332 lines of `OpenBLAS Warning : Detect OpenMP Loop and this
   application may hang` in one 52 s run: a 22 MB `crest.log`. Across about 70
   searches that is gigabytes, and the warning itself mentions hangs. **Fixed (not
   chemistry):** CREST jobs now default to `OPENBLAS_NUM_THREADS=1` unless you
   set it yourself. With that, the log is 33 KB. Alternatively, install the
   OpenMP build: `conda install -c conda-forge "libopenblas=*=*openmp*"`.
3. **`--quick` placed after `--ewin` silently resets the window to 5 kcal/mol**
   (`sorting energy window (EWIN) : 5.0000` in
   `tests/fixtures/crest3/maa_water_quick/crest.quick_after_ewin.log`). Placed
   before `--ewin`, 10 is kept. The pipeline does not use `--quick`. If you add it
   for a test run, put it first.
4. The ensembles vary from run to run. `--quick` kept 4 conformers within 1.1 kcal/mol.
   The full run kept 2. A second, identical full run (inside the SLURM test) kept 3.
   The metadynamics is stochastic, so expect the conformer count, and the number
   of `sp` jobs, to change between reruns. Reruns that matter should restart from
   the same `crest_conformers.xyz` rather than rerunning CREST.

## Places I was unsure of the format

### CREST (genuine output, so these are small)

- Only CREST 3.0.2 (conda-forge build, 2025-08) was checked. A cluster module or the
  static GitHub binary of another 3.0.x should write the same format, but I have not seen one.
- `--nci` mode was not run. I assume it writes the same `crest_conformers.xyz` format.
- These runs used CREST 3's built-in tblite engine (the default; the log says
  `xTB calculation via tblite lib`). `--legacy` mode, which calls the xtb binary,
  was not tried.
- Element symbols are capitalised (`C`, `H`) in CREST 3. ORCA accepts either case.

### ORCA 6 (hand-made fixtures; the bigger risk)

No ORCA output was available. The fixtures follow my memory of ORCA 5/6 and every
one says "HAND-MADE FIXTURE" near the top. Specifically:

1. **`FINAL SINGLE POINT ENERGY`**: confident about the line. I am not sure whether
   ORCA 6 with CPCM prints a second line containing the same phrase (for example a
   gas-phase variant). If it does, `final_energy` takes the last occurrence, which
   might be the wrong one.
2. **Opt+Freq**: I assumed the last `FINAL SINGLE POINT ENERGY` is the "final energy
   evaluation at the stationary point", and that the frequency part prints none.
3. **`G-E(el)                           ...      0.06674219 Eh`**: confident for ORCA 5,
   assumed unchanged in 6. Also unsure whether ORCA 6 applies quasi-RRHO to the entropy
   by default (the README already says to check).
4. **`***imaginary mode***`**: confident for ORCA 4/5, assumed for 6.
5. **`****ORCA TERMINATED NORMALLY****`**: confident.
6. **Unconverged optimizations.** I don't know whether an Opt that hits its maximum
   number of cycles still ends with `TERMINATED NORMALLY`. If it does, `orca_done` calls
   it done, and an unconverged geometry enters the ensemble. Suggestion, not applied:
   for `opt.out` also require `THE OPTIMIZATION HAS CONVERGED` (the HURRAY block; it
   is in the fixture).
7. **SCF failure and truncation text**: approximate. The parsers only depend on the
   termination line being absent, so the exact wording does not matter.
8. **r2SCAN-3c frequencies with CPCM**: I am not sure ORCA 6 has analytic Hessians for
   r2SCAN (a meta-GGA) with CPCM. If not, every `opt` job either fails at the Freq step
   or silently goes numerical, which is much slower. The first real `opt.out` will
   tell; the fixture says "Analytical frequency calculation".
9. **Optimized geometry file**: assumed `opt.xyz`, from ORCA's `<basename>.xyz` convention.
   Confident.
10. **LED output**: no fixture. See item 6 of "Things I think are wrong".

To close these out, run one real sp, opt, hl and cation job on MAA, copy the outputs
over the files in `tests/fixtures/orca6/`, update the expected energies in
`tests/test_parsers.py`, and rerun `pytest`.

## Code changes (all non-chemistry)

- **CREST done test**: requires the normal-termination line, not just the file (above).
- **Job markers**: while a job runs, `<name>.running` holds `local <pid> <host>` or
  `slurm <jobid>_<task> <host>`; `<name>.queued` marks a job submitted to SLURM
  that has not started. Every stage skips live jobs, so a rerun never starts a second
  copy. `opt` waits while any single point of that species is running or queued, and
  `extras` waits for `opt` in the same way. Otherwise a stage could pick conformers
  from a half-finished ensemble.
- **`--slurm`**: one array per stage, split at `max_array_size` (default 1000). ORCA
  tasks get `--ntasks=nprocs` (MPI ranks for `%pal`). CREST tasks get
  `--cpus-per-task=nprocs`. Memory is 1.25 × `maxcore` per CPU. Partition, account,
  wall times, extra `#SBATCH` lines and `module load` lines go in an optional `"slurm"`
  config block (see `SLURM_DEFAULTS` and the README). `all --slurm` is refused, because
  each stage's job list depends on the previous stage's results.
- **`status`**: read-only. Expected job counts use the same selection code as the
  stages. "Not started" for `opt`/`hl`/`extras` reflects only what is known so far,
  so it grows as earlier stages finish.
- **`extras` no longer overwrites `cp/geom.xyz`** once it exists. Before, every rerun
  replaced the fragment labels the README tells you to add by hand.
- **Guards instead of tracebacks**: an `opt.out` without `G-E(el)` is excluded and
  listed (before, `None + float` crashed `extras` and `analyze`). A CREST structure
  without an energy is skipped with a warning. A wrong `orca`/`crest` path gives a
  one-line error.
- **CLI**: `argparse`. Stages, the default config name and the positional config are
  unchanged; `--slurm`, `--no-submit` and `status` are new.

The computing code is unchanged. The `sp`/`opt`/`extras` selection logic moved
into `sp_selection`, `opt_selection` and `best_conformer`, so `status` can reuse it,
with the same arithmetic. `git diff pipeline.py` shows this.

## SLURM caveats (untested on a real cluster)

- ORCA needs its MPI on `PATH`/`LD_LIBRARY_PATH` inside the job. Put the module loads
  in `"slurm": {"setup": [...]}`.
- Jobs run directly in `work/...` on the shared filesystem. ORCA's scratch I/O there can
  be slow. Staging to `$TMPDIR` would be a later improvement.
- The default wall times (`crest` 24 h, `sp` 4 h, `opt` 48 h, `hl` 12 h, `extras` 24 h) are guesses.
- `status` asks `squeue` whether a marked job is still alive. Without `squeue`, a task that
  SLURM killed before it could remove its marker keeps showing as running.
- A job started locally on another machine sharing the filesystem shows as running.
  Its PID cannot be checked from here.
