import csv
from pathlib import Path

FIELDS=("enzyme","source","sequence_position","wild_type","mutant","source_score","source_direction","selected")

def normalize_scores(path: Path, enzyme: str, source: str, score_direction="higher_is_better"):
    with path.open(encoding="utf-8-sig", newline="") as fh: rows=list(csv.DictReader(fh))
    if not rows: raise ValueError(f"empty online result: {path}")
    aliases={"position":"sequence_position","pos":"sequence_position","wt":"wild_type","aa_wt":"wild_type","aa_mut":"mutant","mutation":"mutant","score":"source_score","ev_score":"source_score"}
    normalized=[]
    for raw in rows:
        row={aliases.get(k,k):v for k,v in raw.items()};
        if "sequence_position" not in row or "source_score" not in row: raise ValueError("online result requires position and score columns")
        normalized.append({"enzyme":enzyme,"source":source,"sequence_position":int(row["sequence_position"]),"wild_type":row.get("wild_type", ""),"mutant":row.get("mutant", ""),"source_score":float(row["source_score"]),"source_direction":score_direction,"selected":True})
    try:
        import pandas as pd
        return pd.DataFrame(normalized, columns=FIELDS)
    except ImportError:
        return normalized

