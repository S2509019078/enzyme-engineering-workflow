from pathlib import Path
import csv, re

def parse_average_fxout(path: Path):
    rows=[]
    lines=path.read_text(errors="replace").splitlines()
    for line in lines:
        if not line.strip() or line.lower().startswith("pdb") or line.startswith("#"): continue
        fields=line.split("\t")
        if len(fields)<2: continue
        try: rows.append({"mutation":fields[0].strip(), "ddg_kcal_mol":float(fields[1])})
        except ValueError: continue
    try:
        import pandas as pd
        return pd.DataFrame(rows)
    except ImportError:
        return rows

def classify_ddg(value, stable_max=-1.0, neutral_max=1.0):
    if value <= stable_max: return "stabilizing"
    if value <= neutral_max: return "neutral"
    return "destabilizing"

def build_workbook(frame, path):
    if not hasattr(frame, "to_excel"): raise RuntimeError("pandas and openpyxl are required for XLSX reports")
    path.parent.mkdir(parents=True, exist_ok=True); frame.to_excel(path, index=False)

def build_figures(frame, output_dir):
    if not hasattr(frame, "pivot_table"): raise RuntimeError("pandas, matplotlib and seaborn are required for figures")
    import matplotlib.pyplot as plt, seaborn as sns
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix=frame.pivot_table(index="mutation", values="ddg_kcal_mol")
    ax=sns.heatmap(matrix, cmap="RdBu_r", center=0); ax.figure.savefig(output_dir/"foldx_ddg_heatmap.png", dpi=220, bbox_inches="tight"); plt.close(ax.figure)

