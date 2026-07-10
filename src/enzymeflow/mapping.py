from pathlib import Path
from .models import ResidueMapping

AA3={"ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y"}

def read_fasta(path: Path):
    return "".join(line.strip() for line in path.read_text(encoding="utf-8").splitlines() if not line.startswith(">"))

def pdb_chain_sequence(path: Path, chain: str):
    seen={};
    for line in path.read_text(errors="replace").splitlines():
        if line[:6].strip()!="ATOM" or line[21:22].strip()!=chain: continue
        try: seen[int(line[22:26])]=AA3.get(line[17:20].strip(), "X")
        except ValueError: pass
    return seen

def map_positions(fasta, structure):
    result=[]
    for seq_pos, residue in enumerate(fasta, 1):
        if seq_pos not in structure: continue
        result.append(ResidueMapping(seq_pos,residue,"",seq_pos,residue))
    return result

def validate_candidate_mapping(candidates, mappings):
    mapped={m.sequence_position for m in mappings}
    return [c for c in candidates if c.sequence_position in mapped]

