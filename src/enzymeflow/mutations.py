from .models import Candidate, ResidueMapping

AA = "ACDEFGHIKLMNPQRSTVWY"

def expand_saturation(sites, sequence):
    rows=[]
    for pos in sites:
        wt=sequence[pos]
        if wt not in AA: raise ValueError(f"invalid wild-type residue at {pos}: {wt}")
        rows.extend(Candidate(pos, wt, aa) for aa in AA if aa != wt)
    return rows

def to_foldx_code(candidate, mapping):
    m=mapping.get(candidate.sequence_position)
    if m is None: raise ValueError(f"unmapped sequence position: {candidate.sequence_position}")
    return f"{m.chain}{candidate.mutant}{m.structure_position}{m.chain}"

def write_individual_list(candidates, mapping, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(to_foldx_code(c,mapping) for c in candidates) + "\n", encoding="utf-8")

