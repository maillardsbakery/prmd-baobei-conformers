"""Parsers against genuine CREST 3.0.2 output and hand-made ORCA 6 output.

See tests/fixtures/README.md for where each fixture came from.
"""
import shutil

import pytest

import pipeline
from conftest import FIX

CREST_FULL = FIX / "crest3" / "maa_water_full"
CREST_QUICK = FIX / "crest3" / "maa_water_quick"
ORCA = FIX / "orca6"


# ---------------------------------------------------------------- read_multixyz
def test_read_multixyz_genuine_crest3_conformers():
    ens = pipeline.read_multixyz(CREST_FULL / "crest_conformers.xyz")
    assert [round(e, 8) for e, _ in ens] == [-19.74775542, -19.74703985]
    for e, block in ens:
        lines = block.splitlines()
        assert lines[0] == "12"
        assert len(lines) == 2 + 12
        assert lines[2].split()[0] == "C"
        assert all(len(l.split()) == 4 for l in lines[2:])


def test_read_multixyz_matches_crest_energies_file():
    """Relative energies agree with CREST's own crest.energies (kcal/mol)."""
    for d in (CREST_FULL, CREST_QUICK):
        ens = pipeline.read_multixyz(d / "crest_conformers.xyz")
        rel = [float(l.split()[1]) for l in (d / "crest.energies").read_text().splitlines()]
        assert len(ens) == len(rel)
        e0 = ens[0][0]
        for (e, _), r in zip(ens, rel):
            assert (e - e0) * pipeline.H2KCAL == pytest.approx(r, abs=2e-3)


def test_read_multixyz_sorted_ascending():
    ens = pipeline.read_multixyz(CREST_QUICK / "crest_conformers.xyz")
    energies = [e for e, _ in ens]
    assert len(ens) == 4
    assert energies == sorted(energies)


def test_read_multixyz_best_is_first_conformer():
    best = pipeline.read_multixyz(CREST_FULL / "crest_best.xyz")
    first = pipeline.read_multixyz(CREST_FULL / "crest_conformers.xyz")[0]
    assert len(best) == 1
    assert best[0][0] == first[0]


def test_read_multixyz_rotamer_comment_takes_first_number():
    # crest_rotamers.xyz comments read "  -19.747755309999999  0.1987...  !"
    ens = pipeline.read_multixyz(CREST_QUICK / "crest_rotamers.xyz")
    assert ens[0][0] == pytest.approx(-19.74775531, abs=1e-8)
    assert all(-20 < e < -19 for e, _ in ens)


def test_read_multixyz_block_is_valid_xyz_for_orca(tmp_path):
    _, block = pipeline.read_multixyz(CREST_FULL / "crest_conformers.xyz")[1]
    f = tmp_path / "conf.xyz"
    f.write_text(block)
    again = pipeline.read_multixyz(f)
    assert again == [(pytest.approx(-19.74703985), block)]


def test_read_multixyz_edge_cases(tmp_path):
    f = tmp_path / "x.xyz"
    f.write_text("\n  2\n energy: -1.50000000 gnorm: 0.0001 xtb: 6.7.1\n"
                 "H 0 0 0\nH 0 0 0.74\n\n\n"
                 "2\nno energy here\nH 0 0 0\nH 0 0 0.75")   # no final newline
    ens = pipeline.read_multixyz(f)
    assert ens[0][0] == -1.5
    assert ens[1][0] is None
    assert ens[1][1].endswith("H 0 0 0.75\n")


def test_sp_selection_window():
    cfg = {"sp_window": 0.5}
    sel, total = pipeline.sp_selection(cfg, CREST_QUICK / "crest_conformers.xyz")
    assert total == 4
    assert [k for k, _, _ in sel] == [0, 1]       # 0.000 and 0.451 kcal/mol
    cfg = {"sp_window": 10.0}
    assert len(pipeline.sp_selection(cfg, CREST_QUICK / "crest_conformers.xyz")[0]) == 4


# ---------------------------------------------------------------- final_energy
@pytest.mark.parametrize("name, energy", [
    ("sp_r2scan3c.out", -306.388130517431),
    ("opt_freq_r2scan3c.out", -306.388241226941),     # the last of 4 hits
    ("opt_freq_imaginary.out", -306.388241226941),
    ("hl_wb97mv.out", -306.542187734519),
    ("cation_wb97mv.out", -306.240079918204),
])
def test_final_energy(name, energy):
    assert pipeline.final_energy(ORCA / name) == energy


@pytest.mark.parametrize("name", ["sp_scf_failed.out", "sp_truncated.out"])
def test_final_energy_none_without_result(name):
    assert pipeline.final_energy(ORCA / name) is None


def test_final_energy_takes_last_of_optimization():
    text = (ORCA / "opt_freq_r2scan3c.out").read_text()
    assert text.count("FINAL SINGLE POINT ENERGY") == 4
    assert pipeline.final_energy(ORCA / "opt_freq_r2scan3c.out") != -306.388130517431


# ---------------------------------------------------------------- g_minus_eel
def test_g_minus_eel():
    assert pipeline.g_minus_eel(ORCA / "opt_freq_r2scan3c.out") == 0.06674219


def test_g_minus_eel_consistent_with_final_gibbs():
    out = ORCA / "opt_freq_r2scan3c.out"
    g = float(next(l for l in out.read_text().splitlines()
                   if l.startswith("Final Gibbs free energy")).split()[-2])
    assert pipeline.final_energy(out) + pipeline.g_minus_eel(out) == pytest.approx(g, abs=1e-8)


@pytest.mark.parametrize("name", ["sp_r2scan3c.out", "hl_wb97mv.out", "sp_truncated.out"])
def test_g_minus_eel_none_without_freq(name):
    assert pipeline.g_minus_eel(ORCA / name) is None


# ---------------------------------------------------------------- has_imaginary
@pytest.mark.parametrize("name, expected", [
    ("opt_freq_r2scan3c.out", False),
    ("opt_freq_imaginary.out", True),
    ("sp_r2scan3c.out", False),
])
def test_has_imaginary(name, expected):
    assert pipeline.has_imaginary(ORCA / name) is expected


# ---------------------------------------------------------------- orca_done
@pytest.mark.parametrize("name, expected", [
    ("sp_r2scan3c.out", True),
    ("opt_freq_r2scan3c.out", True),
    ("opt_freq_imaginary.out", True),
    ("hl_wb97mv.out", True),
    ("cation_wb97mv.out", True),
    ("sp_scf_failed.out", False),
    ("sp_truncated.out", False),
])
def test_orca_done(name, expected):
    assert pipeline.orca_done(ORCA / name) is expected


def test_orca_done_missing_file(tmp_path):
    assert pipeline.orca_done(tmp_path / "nope.out") is False


# ---------------------------------------------------------------- oxidation_potential
def _ox_tree(root, hl_files, cation_file):
    base = root / "work" / "MAA" / "neutral"
    for k, f in enumerate(hl_files):
        d = base / "hl" / f"c{k:03d}"
        d.mkdir(parents=True)
        shutil.copy(ORCA / f, d / "hl.out")
    if cation_file:
        d = base / "extras" / "cation"
        d.mkdir(parents=True)
        shutil.copy(ORCA / cation_file, d / "cation.out")


OX_CFG = {"species": {"MAA": {"variants": {"neutral": {}}}}}


def test_oxidation_potential(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _ox_tree(tmp_path, ["hl_wb97mv.out"], "cation_wb97mv.out")
    expected = (-306.240079918204 - -306.542187734519) * 27.2114 - 4.28
    assert pipeline.oxidation_potential(OX_CFG, "MAA") == pytest.approx(expected)
    assert expected == pytest.approx(3.9408, abs=1e-4)


def test_oxidation_potential_she_reference(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _ox_tree(tmp_path, ["hl_wb97mv.out"], "cation_wb97mv.out")
    a = pipeline.oxidation_potential(OX_CFG, "MAA")
    b = pipeline.oxidation_potential({**OX_CFG, "she_ref_V": 4.44}, "MAA")
    assert a - b == pytest.approx(0.16)


def test_oxidation_potential_uses_lowest_neutral(tmp_path, monkeypatch):
    # c000 is the lower wB97M-V energy; c001 (r2SCAN-3c, higher) must not be used
    monkeypatch.chdir(tmp_path)
    _ox_tree(tmp_path, ["hl_wb97mv.out", "sp_r2scan3c.out"], "cation_wb97mv.out")
    expected = (-306.240079918204 - -306.542187734519) * 27.2114 - 4.28
    assert pipeline.oxidation_potential(OX_CFG, "MAA") == pytest.approx(expected)


def test_oxidation_potential_skips_failed_neutral(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _ox_tree(tmp_path, ["sp_truncated.out", "hl_wb97mv.out"], "cation_wb97mv.out")
    assert pipeline.oxidation_potential(OX_CFG, "MAA") == pytest.approx(3.9408, abs=1e-4)


@pytest.mark.parametrize("hl, cation", [
    (["hl_wb97mv.out"], None),
    (["hl_wb97mv.out"], "sp_scf_failed.out"),
    ([], "cation_wb97mv.out"),
    (["sp_truncated.out"], "cation_wb97mv.out"),
])
def test_oxidation_potential_none_when_incomplete(tmp_path, monkeypatch, hl, cation):
    monkeypatch.chdir(tmp_path)
    _ox_tree(tmp_path, hl, cation)
    assert pipeline.oxidation_potential(OX_CFG, "MAA") is None


def test_oxidation_potential_unknown_species(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert pipeline.oxidation_potential(OX_CFG, "nope") is None


# ---------------------------------------------------------------- crest_done
def test_crest_done_genuine_finished_run():
    assert pipeline.crest_done(CREST_FULL)
    assert pipeline.crest_done(CREST_QUICK)


def test_crest_done_false_mid_run():
    # snapshot taken while CREST 3.0.2 was still running: the conformer file
    # already exists (3 structures; the finished run kept 2)
    mid = FIX / "crest3" / "maa_water_midrun"
    assert len(pipeline.read_multixyz(mid / "crest_conformers.xyz")) == 3
    assert not pipeline.crest_done(mid)


def test_crest_done_needs_log(tmp_path):
    shutil.copy(CREST_FULL / "crest_conformers.xyz", tmp_path)
    assert not pipeline.crest_done(tmp_path)
    shutil.copy(CREST_FULL / "crest.log", tmp_path)
    assert pipeline.crest_done(tmp_path)
