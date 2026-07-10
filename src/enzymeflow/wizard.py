from __future__ import annotations

from datetime import datetime
from pathlib import Path
import csv
import re
import shutil
import yaml


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    return cleaned or "enzyme_run"


def create_run_directory(base: Path, project_name: str, now: datetime | None = None) -> Path:
    stamp = (now or datetime.now()).strftime("%Y%m%d_%H%M%S")
    run_dir = base / f"{stamp}_{safe_name(project_name)}"
    suffix = 1
    while run_dir.exists():
        run_dir = base / f"{stamp}_{safe_name(project_name)}_{suffix:02d}"
        suffix += 1
    for relative in ("config", "inputs/fasta", "inputs/pdb", "inputs/online", "work", "results", "logs"):
        (run_dir / relative).mkdir(parents=True, exist_ok=True)
    return run_dir


def detect_pdb_chains(pdb_path: Path) -> list[str]:
    chains = []
    seen = set()
    for line in pdb_path.read_text(errors="replace").splitlines():
        if line[:6].strip() != "ATOM":
            continue
        chain = line[21:22].strip()
        if chain not in seen:
            seen.add(chain)
            chains.append(chain)
    return chains


def copy_input(source: Path, destination: Path) -> Path:
    source = source.expanduser()
    if not source.is_file():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return destination


def resolve_foldx(value: str) -> str:
    text = value.strip() or "foldx"
    explicit = Path(text).expanduser()
    if explicit.is_file():
        return str(explicit.resolve())
    found = shutil.which(text)
    if found:
        return str(Path(found).resolve())
    return text


def write_project(run_dir: Path, row: dict, foldx_path: str = "foldx") -> Path:
    enzymes_path = run_dir / "config" / "enzymes.tsv"
    columns = ["name", "fasta", "pdb", "chain", "online_result", "online_mode", "online_source", "score_direction", "selection_threshold", "saturation_scan", "batch_size", "number_of_runs", "stable_max_ddg", "neutral_max_ddg"]
    with enzymes_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerow({key: row.get(key, "") for key in columns})
    config = {
        "tools": {"foldx": foldx_path},
        "paths": {"enzymes": "config/enzymes.tsv", "work": "work", "results": "results"},
        "settings": {
            "online_mode": "mutation_effect",
            "online_source": "online_model",
            "score_direction": "higher_is_better",
            "selection_threshold": 0.0,
            "saturation_scan": True,
            "number_of_runs": 1,
            "out_pdb": False,
            "stable_max_ddg": -1.0,
            "neutral_max_ddg": 1.0,
            "batch_size": 200,
        },
    }
    config_path = run_dir / "config" / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False, allow_unicode=True), encoding="utf-8")
    (run_dir / "RUN_INFO.txt").write_text("本目录是一次独立酶改造运行。输入、配置、中间文件、FoldX结果和日志不会与其他运行混用。\n", encoding="utf-8")
    return config_path


def interactive_wizard(input_fn=input, print_fn=print, base_dir: Path = Path("runs")) -> Path:
    print_fn("EnzymeFlow 交互式项目向导")
    project_name = input_fn("项目名称: ").strip() or "enzyme"
    run_dir = create_run_directory(base_dir, project_name)

    fasta_source = Path(input_fn("FASTA 文件路径: ").strip())
    pdb_source = Path(input_fn("PDB 文件路径: ").strip())
    online_source = Path(input_fn("在线预测结果 CSV/TSV 路径: ").strip())

    fasta_dest = copy_input(fasta_source, run_dir / "inputs" / "fasta" / fasta_source.name)
    pdb_dest = copy_input(pdb_source, run_dir / "inputs" / "pdb" / pdb_source.name)
    online_dest = copy_input(online_source, run_dir / "inputs" / "online" / online_source.name)

    chains = detect_pdb_chains(pdb_dest)
    print_fn("检测到 PDB 链: " + (", ".join(chain or "<空链>" for chain in chains) if chains else "未检测到"))
    default_chain = chains[0] if len(chains) == 1 else ""
    chain = input_fn(f"选择用于映射的链{f' [{default_chain}]' if default_chain else ''}: ").strip() or default_chain
    if chains and chain not in chains:
        raise ValueError(f"链 {chain!r} 不在检测到的链中: {chains}")

    direction = input_fn("在线分数方向 higher_is_better/lower_is_better [higher_is_better]: ").strip() or "higher_is_better"
    if direction not in {"higher_is_better", "lower_is_better"}:
        raise ValueError("score_direction 必须是 higher_is_better 或 lower_is_better")
    threshold_text = input_fn("筛选阈值 [0.0]: ").strip()
    threshold = float(threshold_text) if threshold_text else 0.0
    saturation = input_fn("对入选位点做19种饱和突变？ Y/n [Y]: ").strip().lower() not in {"n", "no", "0"}
    foldx_path = resolve_foldx(input_fn("FoldX 可执行文件路径或命令名 [foldx]: "))

    row = {
        "name": safe_name(project_name),
        "fasta": str(fasta_dest.relative_to(run_dir)),
        "pdb": str(pdb_dest.relative_to(run_dir)),
        "chain": chain,
        "online_result": str(online_dest.relative_to(run_dir)),
        "online_mode": "mutation_effect",
        "online_source": "online_model",
        "score_direction": direction,
        "selection_threshold": threshold,
        "saturation_scan": str(saturation).lower(),
        "batch_size": 200,
        "number_of_runs": 1,
        "stable_max_ddg": -1.0,
        "neutral_max_ddg": 1.0,
    }
    config_path = write_project(run_dir, row, foldx_path=foldx_path)
    print_fn(f"项目已创建: {run_dir}")
    print_fn(f"先检查: enzymeflow check --config {config_path}")
    print_fn(f"运行: enzymeflow all --config {config_path}")
    return config_path
