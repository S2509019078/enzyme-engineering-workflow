import csv
import re
from pathlib import Path

FIELDS = ("enzyme", "source", "sequence_position", "wild_type", "mutant", "source_score", "source_direction", "selected")
MUTATION_RE = re.compile(r"^([A-Za-z])([0-9]+)([A-Za-z])$")


def _pick(row, *names):
    lowered = {str(key).strip().lower(): value for key, value in row.items()}
    for name in names:
        if name in lowered and str(lowered[name]).strip() != "":
            return str(lowered[name]).strip()
    return ""


def normalize_scores(path: Path, enzyme: str, source: str, score_direction="higher_is_better", selection_threshold=0.0):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"empty online result: {path}")

    normalized = []
    for raw in rows:
        mutation = _pick(raw, "mutation", "mutant")
        wild_type = _pick(raw, "wild_type", "wt", "aa_wt").upper()
        mutant = _pick(raw, "mutant_aa", "aa_mut", "alt").upper()
        position_text = _pick(raw, "sequence_position", "position", "pos", "site")
        match = MUTATION_RE.match(mutation.upper()) if mutation else None
        if match:
            wild_type = wild_type or match.group(1)
            position_text = position_text or match.group(2)
            mutant = mutant or match.group(3)
        score_text = _pick(raw, "source_score", "score", "ev_score", "evolutionary_index", "esm_score")
        if not position_text or not score_text:
            raise ValueError("mutation-effect input requires a position/mutation column and a score column")
        score = float(score_text)
        selected = score >= selection_threshold if score_direction == "higher_is_better" else score <= selection_threshold
        normalized.append({"enzyme": enzyme, "source": source, "sequence_position": int(position_text), "wild_type": wild_type, "mutant": mutant, "source_score": score, "source_direction": score_direction, "selected": selected})

    try:
        import pandas as pd
        return pd.DataFrame(normalized, columns=FIELDS)
    except ImportError:
        return normalized


def normalize_pairwise_couplings(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    normalized = []
    for row in rows:
        i = _pick(row, "i", "position_i")
        j = _pick(row, "j", "position_j")
        score = _pick(row, "cn", "score", "coupling_score")
        if not i or not j or not score:
            raise ValueError("pairwise coupling input requires i, j and cn/score columns")
        normalized.append({"position_i": int(i), "position_j": int(j), "coupling_score": float(score)})
    try:
        import pandas as pd
        return pd.DataFrame(normalized)
    except ImportError:
        return normalized
