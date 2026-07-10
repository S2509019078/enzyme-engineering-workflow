import pytest

from enzymeflow.report import classify_ddg, parse_dif_fxout


def test_parse_foldx_dif_and_classify(tmp_path):
    path = tmp_path / "Dif_x.fxout"
    path.write_text(
        "Pdb\tTotal Energy\tBackbone Hbond\n"
        "x_1\t-2.5\t0.1\n"
        "x_2\t1.2\t0.2\n",
        encoding="utf-8",
    )
    frame = parse_dif_fxout(path, ["AA1V;", "AA1G;"])
    assert list(frame["mutation"]) == ["AA1V;", "AA1G;"]
    assert list(frame["ddg_kcal_mol"]) == [-2.5, 1.2]
    assert list(frame["mutant"]) == ["V", "G"]
    assert classify_ddg(-2.5, -1.0, 1.0) == "stabilizing"


def test_parse_dif_requires_total_energy(tmp_path):
    path = tmp_path / "Dif_x.fxout"
    path.write_text("Pdb\tOther\nx_1\t-2.5\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Total Energy"):
        parse_dif_fxout(path)


def test_mutation_count_must_match_foldx_rows(tmp_path):
    path = tmp_path / "Dif_x.fxout"
    path.write_text("Pdb\tTotal Energy\nx_1\t-2.5\nx_2\t1.2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="mutation count"):
        parse_dif_fxout(path, ["AA1V;"])
