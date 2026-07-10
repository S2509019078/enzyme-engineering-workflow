from pathlib import Path
import csv

def load_enzymes(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as fh: rows=list(csv.DictReader(fh, delimiter="\t"))
    required={"name","fasta","pdb","chain"}
    if not rows or required-set(rows[0]): raise ValueError(f"enzyme table requires: {', '.join(sorted(required))}")
    names=set(); result=[]
    for row in rows:
        if row["name"] in names: raise ValueError(f"duplicate enzyme: {row['name']}")
        names.add(row["name"]); result.append({k:(v or "").strip() for k,v in row.items()})
    return result

