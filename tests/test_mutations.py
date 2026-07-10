from enzymeflow.mutations import expand_saturation, to_foldx_code
from enzymeflow.models import ResidueMapping

def test_each_site_gets_nineteen_non_wildtype_mutations():
    rows = expand_saturation([1], {1: "A"})
    assert len(rows) == 19
    assert all(r.mutant != "A" for r in rows)

def test_foldx_code_uses_structure_residue_number():
    c = expand_saturation([1], {1: "A"})[0]
    mapping = {1: ResidueMapping(1, "A", "A", 429, "A")}
    assert to_foldx_code(c, mapping).startswith("AA429")

