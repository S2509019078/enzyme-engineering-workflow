import pytest

from enzymeflow.models import ResidueMapping
from enzymeflow.mutations import expand_saturation, partition_candidates, to_foldx_code


def test_each_site_gets_nineteen_non_wildtype_mutations():
    rows = expand_saturation([1], "A")
    assert len(rows) == 19
    assert all(row.mutant != "A" for row in rows)


def test_expand_saturation_uses_one_based_positions():
    rows = expand_saturation([2], "MA")
    assert all(row.wild_type == "A" for row in rows)


def test_foldx_code_uses_wildtype_chain_structure_number_and_mutant():
    candidate = next(row for row in expand_saturation([1], "A") if row.mutant == "V")
    mapping = {1: ResidueMapping(1, "A", "A", 429, "A")}
    assert to_foldx_code(candidate, mapping) == "AA429V;"


def test_foldx_code_rejects_candidate_fasta_mismatch():
    candidate = next(row for row in expand_saturation([1], "A") if row.mutant == "V")
    mapping = {1: ResidueMapping(1, "G", "A", 429, "G")}
    with pytest.raises(ValueError, match="candidate_fasta_mismatch"):
        to_foldx_code(candidate, mapping)


def test_foldx_code_rejects_fasta_pdb_mismatch():
    candidate = next(row for row in expand_saturation([1], "A") if row.mutant == "V")
    mapping = {1: ResidueMapping(1, "A", "A", 429, "G")}
    with pytest.raises(ValueError, match="fasta_pdb_mismatch"):
        to_foldx_code(candidate, mapping)


def test_partition_candidates_reports_rejection_reason():
    candidate = next(row for row in expand_saturation([1], "A") if row.mutant == "V")
    mapping = {1: ResidueMapping(1, "A", "A", 429, "G")}
    valid, rejected = partition_candidates([candidate], mapping)
    assert valid == []
    assert rejected[0][1] == "fasta_pdb_mismatch"
