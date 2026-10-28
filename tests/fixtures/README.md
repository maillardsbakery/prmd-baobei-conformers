# Test fixtures

## crest3/: genuine CREST 3.0.2 output

These files came from real runs in the cloud session that wrote these tests:
CREST 3.0.2 (commit 6914a25) and xtb 6.7.1 from conda-forge, on methacrylic acid
(MAA) in GFN2-xTB/ALPB(water). The start geometry (`maa_water_full/start.xyz`)
was built by hand and pre-optimized with `xtb --opt --gfn 2 --alpb water`.

| Directory | Command | Notes |
| --- | --- | --- |
| `maa_water_full/` | `crest start.xyz --gfn2 --alpb water --chrg 0 --uhf 0 --T 4 --ewin 10.0` | Run by `pipeline.py crest`, so this is exactly the pipeline's command line. 2 conformers, 76 s on 4 threads. |
| `maa_water_quick/` | same with `--quick` placed **before** `--ewin` | 4 conformers. `crest.quick_after_ewin.log` is an earlier run with `--quick` placed **after** `--ewin`: there CREST silently used a 5 kcal/mol window (see "sorting energy window"). |
| `maa_water_midrun/` | the `maa_water_full` command, copied while CREST was still running | CREST 3 rewrites `crest_conformers.xyz` after every cycle, so the file exists before the run ends. This snapshot has 3 structures; the finished run had 2. |

`crest.quick_after_ewin.log` had 184,332 lines of
`OpenBLAS Warning : Detect OpenMP Loop and this application may hang...`
removed (22 MB). The warnings come from conda-forge's pthreads OpenBLAS
build and stop with `OPENBLAS_NUM_THREADS=1`, which the other runs used.
Nothing else in these files was edited.

## orca6/: hand-made ORCA 6 output

ORCA could not be installed here (it needs a licensed download), so these are
**written by hand** from memory of ORCA 5/6 output. Every file says so near its
top. Sections pipeline.py never reads are shortened or left out. The lines the
parsers depend on follow the real format as closely as I could get them:

- `FINAL SINGLE POINT ENERGY      -306.388130517431`
- `G-E(el)                           ...      0.06674219 Eh     41.88 kcal/mol`
- `     6:       -38.17 cm**-1 ***imaginary mode***`
- `                             ****ORCA TERMINATED NORMALLY****`

REVIEW.md lists what I was unsure of. Before trusting the pipeline on a cluster,
run one real job of each kind (sp, opt, hl, cation), copy the outputs over these
files, and rerun `pytest`. The energies in the tests will need updating; the
parsers should not.

| File | Job |
| --- | --- |
| `sp_r2scan3c.out` | r2SCAN-3c CPCM(water) single point (sp stage) |
| `opt_freq_r2scan3c.out` | r2SCAN-3c Opt Freq, 3 cycles + final evaluation, no imaginary modes (opt stage) |
| `opt_freq_imaginary.out` | same, with one imaginary mode |
| `hl_wb97mv.out` | wB97M-V/def2-TZVP single point (hl stage) |
| `cation_wb97mv.out` | wB97M-V UKS doublet, vertical radical cation (extras stage) |
| `sp_scf_failed.out` | SCF not converged, error termination |
| `sp_truncated.out` | killed mid-SCF (wall-time or OOM): no termination line |

The energies are made up but self-consistent: the SP total = SCF energy + D4 + gCP,
and in the Opt+Freq files E(el) + G-E(el) = the final Gibbs free energy.

`tests/fake_bin/crest` and `tests/fake_bin/orca` emit these files, with the job's
input and a made-up energy spliced in.
