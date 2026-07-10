from enzymeflow.mapping import map_positions


def test_alignment_maps_sequence_positions_to_shifted_structure_numbers():
    fasta = "MQQACDEFG"
    structure = [(213, "", "A"), (214, "", "C"), (215, "", "D"), (216, "", "E"), (217, "", "F"), (218, "", "G")]
    mapping = map_positions(fasta, structure, chain="A")
    by_sequence = {item.sequence_position: item for item in mapping}
    assert by_sequence[4].structure_position == 213
    assert by_sequence[4].chain == "A"
    assert by_sequence[9].structure_position == 218


def test_alignment_skips_unresolved_structure_residues_without_assuming_same_number():
    fasta = "ACDEFG"
    structure = [(10, "", "A"), (11, "", "C"), (13, "", "E"), (14, "", "F"), (15, "", "G")]
    mapping = map_positions(fasta, structure, chain="B")
    pairs = {(item.sequence_position, item.structure_position) for item in mapping}
    assert (1, 10) in pairs
    assert (2, 11) in pairs
    assert (4, 13) in pairs
    assert (6, 15) in pairs
