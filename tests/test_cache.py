from enzymeflow.cache import build_payload, manifest_valid, write_manifest


def test_manifest_invalidates_when_mutation_list_changes(tmp_path):
    individual = tmp_path / "individual_list.txt"
    repaired = tmp_path / "enzyme_Repair.pdb"
    output = tmp_path / "Dif_enzyme_Repair.fxout"
    manifest = tmp_path / "manifest.json"
    individual.write_text("AA1V;\n", encoding="utf-8")
    repaired.write_text("ATOM\n", encoding="utf-8")
    output.write_text("Pdb\tTotal Energy\nx\t-1.0\n", encoding="utf-8")
    payload = build_payload({"individual_list": individual, "repaired_pdb": repaired}, {"number_of_runs": 1})
    write_manifest(manifest, payload, [output])
    assert manifest_valid(manifest, payload, [output])
    individual.write_text("AA1G;\n", encoding="utf-8")
    changed = build_payload({"individual_list": individual, "repaired_pdb": repaired}, {"number_of_runs": 1})
    assert not manifest_valid(manifest, changed, [output])


def test_manifest_rejects_missing_output(tmp_path):
    source = tmp_path / "source.txt"
    output = tmp_path / "output.txt"
    manifest = tmp_path / "manifest.json"
    source.write_text("input", encoding="utf-8")
    output.write_text("result", encoding="utf-8")
    payload = build_payload({"source": source})
    write_manifest(manifest, payload, [output])
    output.unlink()
    assert not manifest_valid(manifest, payload, [output])
