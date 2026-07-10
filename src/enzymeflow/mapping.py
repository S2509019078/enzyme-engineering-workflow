from pathlib import Path

from Bio.Align import PairwiseAligner

from .models import ResidueMapping

AA3 = {"ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y","MSE":"M"}


def read_fasta(path: Path):
    sequence = "".join(line.strip() for line in path.read_text(encoding="utf-8").splitlines() if not line.startswith(">"))
    sequence = sequence.replace(" ", "").upper()
    if not sequence:
        raise ValueError(f"empty FASTA: {path}")
    return sequence


def pdb_chain_residues(path: Path, chain: str):
    residues = []
    seen = set()
    for line in path.read_text(errors="replace").splitlines():
        if line[:6].strip() != "ATOM" or line[21:22].strip() != chain:
            continue
        try:
            position = int(line[22:26])
        except ValueError:
            continue
        insertion_code = line[26:27].strip()
        key = (position, insertion_code)
        if key in seen:
            continue
        seen.add(key)
        residues.append((position, insertion_code, AA3.get(line[17:20].strip().upper(), "X")))
    if not residues:
        raise ValueError(f"no residues found for chain {chain!r} in {path}")
    return residues


def map_positions(fasta: str, structure_residues, chain: str = ""):
    structure_sequence = "".join(item[2] for item in structure_residues)
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(fasta, structure_sequence)[0]

    result = []
    for (seq_start, seq_end), (str_start, str_end) in zip(alignment.aligned[0], alignment.aligned[1]):
        length = min(seq_end - seq_start, str_end - str_start)
        for offset in range(length):
            sequence_position = seq_start + offset + 1
            structure_index = str_start + offset
            structure_position, _icode, structure_residue = structure_residues[structure_index]
            result.append(ResidueMapping(sequence_position, fasta[sequence_position - 1], chain, structure_position, structure_residue))
    return result


def build_mapping(fasta_path: Path, pdb_path: Path, chain: str):
    fasta = read_fasta(fasta_path)
    return map_positions(fasta, pdb_chain_residues(pdb_path, chain), chain=chain)


def validate_candidate_mapping(candidates, mappings):
    mapping = {item.sequence_position: item for item in mappings}
    valid = []
    for candidate in candidates:
        item = mapping.get(candidate.sequence_position)
        if item and item.sequence_residue == candidate.wild_type and item.structure_residue == candidate.wild_type:
            valid.append(candidate)
    return valid
