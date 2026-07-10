from pathlib import Path

from .models import Candidate

AA = "ACDEFGHIKLMNPQRSTVWY"


def expand_saturation(sites, sequence):
    rows = []
    for pos in sites:
        if pos < 1 or pos > len(sequence):
            raise ValueError(f"sequence position out of range: {pos}")
        wt = sequence[pos - 1].upper()
        if wt not in AA:
            raise ValueError(f"invalid wild-type residue at {pos}: {wt}")
        rows.extend(Candidate(pos, wt, aa) for aa in AA if aa != wt)
    return rows


def to_foldx_code(candidate, mapping):
    item = mapping.get(candidate.sequence_position)
    if item is None:
        raise ValueError(f"unmapped sequence position: {candidate.sequence_position}")
    if item.sequence_residue != candidate.wild_type:
        raise ValueError(f"wild-type mismatch at {candidate.sequence_position}: candidate={candidate.wild_type}, mapped={item.sequence_residue}")
    return f"{item.structure_residue}{item.chain}{item.structure_position}{candidate.mutant};"


def write_individual_list(candidates, mapping, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [to_foldx_code(candidate, mapping) for candidate in candidates]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
