from pathlib import Path
import re

MUTATION_RE = re.compile(r"^([A-Za-z])([A-Za-z])(-?[0-9]+)([A-Za-z]);?$")


def parse_dif_fxout(path: Path, mutations=None):
    lines = [line for line in path.read_text(errors="replace").splitlines() if line.strip() and not line.startswith("#")]
    header = None
    rows = []
    for line in lines:
        fields = [field.strip() for field in line.split("\t")]
        if header is None and any(field.lower() in {"pdb", "total energy"} for field in fields):
            header = fields
            continue
        if header is None or len(fields) < 2:
            continue
        record = dict(zip(header, fields))
        energy_key = next((key for key in header if key.lower() == "total energy"), None)
        if energy_key is None:
            raise ValueError("Dif fxout does not contain a Total Energy column")
        try:
            ddg = float(record[energy_key])
        except ValueError:
            continue
        rows.append({"foldx_row": len(rows) + 1, "pdb": record.get("Pdb", record.get("pdb", "")), "ddg_kcal_mol": ddg})
    if not rows:
        raise ValueError(f"no FoldX difference rows parsed from {path}")
    if mutations is not None:
        mutations = list(mutations)
        if len(mutations) != len(rows):
            raise ValueError(f"mutation count ({len(mutations)}) does not match FoldX rows ({len(rows)})")
        for row, mutation in zip(rows, mutations):
            code = mutation.strip()
            row["mutation"] = code
            match = MUTATION_RE.match(code)
            if match:
                row["wild_type"] = match.group(1).upper()
                row["chain"] = match.group(2)
                row["structure_position"] = int(match.group(3))
                row["mutant"] = match.group(4).upper()
    try:
        import pandas as pd
        return pd.DataFrame(rows)
    except ImportError:
        return rows


def classify_ddg(value, stable_max=-1.0, neutral_max=1.0):
    if value <= stable_max:
        return "stabilizing"
    if value <= neutral_max:
        return "neutral"
    return "destabilizing"


def build_workbook(frame, path):
    if not hasattr(frame, "to_excel"):
        raise RuntimeError("pandas and openpyxl are required for XLSX reports")
    output = frame.copy()
    output["stability_class"] = output["ddg_kcal_mol"].map(classify_ddg)
    path.parent.mkdir(parents=True, exist_ok=True)
    output.to_excel(path, index=False)
    return path


def build_figures(frame, output_dir):
    if not hasattr(frame, "pivot_table"):
        raise RuntimeError("pandas and matplotlib are required for figures")
    import matplotlib.pyplot as plt
    output_dir.mkdir(parents=True, exist_ok=True)
    required = {"mutant", "sequence_position", "ddg_kcal_mol"}
    if not required.issubset(frame.columns):
        raise ValueError(f"heatmap requires columns: {', '.join(sorted(required))}")
    matrix = frame.pivot_table(index="mutant", columns="sequence_position", values="ddg_kcal_mol", aggfunc="mean")
    figure, axis = plt.subplots(figsize=(max(7, 0.5 * len(matrix.columns)), 7))
    image = axis.imshow(matrix.values, aspect="auto", interpolation="nearest")
    axis.set_xticks(range(len(matrix.columns)), labels=[str(value) for value in matrix.columns], rotation=90)
    axis.set_yticks(range(len(matrix.index)), labels=list(matrix.index))
    axis.set_xlabel("Sequence position")
    axis.set_ylabel("Mutant amino acid")
    figure.colorbar(image, ax=axis, label="FoldX ΔΔG (kcal/mol)")
    figure.tight_layout()
    output = output_dir / "foldx_ddg_heatmap.png"
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)
    return output
