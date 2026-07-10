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


def validate_candidate(candidate, mapping):
    item = mapping.get(candidate.sequence_position)
    if item is None:
        return False, "unmapped_sequence_position"
    if candidate.wild_type not in AA or candidate.mutant not in AA:
        return False, "non_standard_amino_acid"
    if candidate.wild_type == candidate.mutant:
        return False, "synonymous_substitution"
    if item.sequence_residue != candidate.wild_type:
        return False, "candidate_fasta_mismatch"
    if item.structure_residue != candidate.wild_type:
        return False, "fasta_pdb_mismatch"
    return True, ""


def to_foldx_code(candidate, mapping):
    valid, reason = validate_candidate(candidate, mapping)
    if not valid:
        raise ValueError(f"invalid candidate at {candidate.sequence_position}: {reason}")
    item = mapping[candidate.sequence_position]
    return f"{item.structure_residue}{item.chain}{item.structure_position}{candidate.mutant};"


def partition_candidates(candidates, mapping):
    valid = []
    rejected = []
    for candidate in candidates:
        ok, reason = validate_candidate(candidate, mapping)
        if ok:
            valid.append(candidate)
        else:
            rejected.append((candidate, reason))
    return valid, rejected


def write_individual_list(candidates, mapping, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [to_foldx_code(candidate, mapping) for candidate in candidates]
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return path
