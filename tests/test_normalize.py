import pytest

from enzymeflow.normalize import normalize_pairwise_couplings, normalize_scores


def test_normalize_mutation_string_and_threshold(tmp_path):
    path = tmp_path / "scores.csv"
    path.write_text("mutation,score\nA12V,1.5\nG20D,-0.2\n", encoding="utf-8")
    frame = normalize_scores(path, "enzyme", "EVmutation", selection_threshold=0.0)
    assert list(frame["sequence_position"]) == [12, 20]
    assert list(frame["wild_type"]) == ["A", "G"]
    assert list(frame["mutant"]) == ["V", "D"]
    assert list(frame["selected"]) == [True, False]


def test_lower_is_better_selection(tmp_path):
    path = tmp_path / "scores.csv"
    path.write_text("position,wt,aa_mut,score\n12,A,V,-1.5\n20,G,D,0.2\n", encoding="utf-8")
    frame = normalize_scores(path, "enzyme", "model", score_direction="lower_is_better", selection_threshold=0.0)
    assert list(frame["selected"]) == [True, False]


def test_pairwise_schema_is_separate(tmp_path):
    path = tmp_path / "couplings.csv"
    path.write_text("i,j,cn\n12,20,0.8\n", encoding="utf-8")
    frame = normalize_pairwise_couplings(path)
    assert frame.iloc[0].to_dict() == {"position_i": 12.0, "position_j": 20.0, "coupling_score": 0.8}


def test_mutation_effect_requires_score(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("mutation\nA12V\n", encoding="utf-8")
    with pytest.raises(ValueError, match="score"):
        normalize_scores(path, "enzyme", "model")
