from dataclasses import dataclass

@dataclass(frozen=True)
class ResidueMapping:
    sequence_position: int; sequence_residue: str; chain: str; structure_position: int; structure_residue: str

@dataclass(frozen=True)
class Candidate:
    sequence_position: int; wild_type: str; mutant: str

@dataclass(frozen=True)
class FoldXResult:
    mutation: str; ddg_kcal_mol: float; status: str = "complete"

