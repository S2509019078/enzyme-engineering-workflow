from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import csv
import shutil
import yaml

import pandas as pd

from .foldx import FoldXRunner
from .mapping import build_mapping, read_fasta
from .models import Candidate
from .mutations import expand_saturation, write_individual_list
from .normalize import normalize_pairwise_couplings, normalize_scores
from .report import build_figures, build_workbook, parse_dif_fxout


class WorkflowConfig:
    def __init__(self, path: Path):
        self.path = path.resolve()
        self.root = self.path.parent.parent
        self.data = yaml.safe_load(self.path.read_text(encoding="utf-8")) or {}
        self.tools = self.data.get("tools", {})
        self.settings = self.data.get("settings", {})
        self.paths = self.data.get("paths", {})

    def resolve(self, value: str | Path) -> Path:
        path = Path(value).expanduser()
        return path if path.is_absolute() else self.root / path

    @property
    def enzymes_table(self):
        return self.resolve(self.paths.get("enzymes", "config/enzymes.tsv"))

    @property
    def work_dir(self):
        return self.resolve(self.paths.get("work", "work"))

    @property
    def results_dir(self):
        return self.resolve(self.paths.get("results", "results"))


class EnzymeWorkflow:
    def __init__(self, config: WorkflowConfig):
        self.config = config
        self.enzymes = self._load_enzymes(config.enzymes_table)

    @staticmethod
    def _load_enzymes(path: Path):
        if not path.exists():
            raise FileNotFoundError(path)
        with path.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        required = {"name", "fasta", "pdb", "chain", "online_result"}
        if not rows:
            raise ValueError(f"empty enzyme table: {path}")
        missing = required - set(rows[0])
        if missing:
            raise ValueError(f"missing enzyme columns: {', '.join(sorted(missing))}")
        return rows

    def selected(self, name=None):
        rows = self.enzymes if not name else [row for row in self.enzymes if row["name"] == name]
        if name and not rows:
            raise ValueError(f"unknown enzyme: {name}")
        return rows

    def check(self):
        problems = []
        if not self.config.enzymes_table.exists():
            problems.append(f"missing enzyme table: {self.config.enzymes_table}")
        foldx = self.config.tools.get("foldx")
        if not foldx or not self.config.resolve(foldx).exists():
            problems.append("FoldX executable not found; configure tools.foldx")
        for row in self.enzymes:
            for key in ("fasta", "pdb", "online_result"):
                path = self.config.resolve(row[key])
                if not path.exists():
                    problems.append(f"missing {key} for {row['name']}: {path}")
        return problems

    def normalize(self, name=None, force=False):
        outputs = []
        threshold = float(self.config.settings.get("selection_threshold", 0.0))
        direction = self.config.settings.get("score_direction", "higher_is_better")
        source = self.config.settings.get("online_source", "mutation_effect")
        mode = self.config.settings.get("online_mode", "mutation_effect")
        for row in self.selected(name):
            out = self.config.work_dir / "normalized" / f"{row['name']}.csv"
            if out.exists() and not force:
                outputs.append(out)
                continue
            out.parent.mkdir(parents=True, exist_ok=True)
            source_path = self.config.resolve(row["online_result"])
            if mode == "pairwise_coupling":
                frame = normalize_pairwise_couplings(source_path)
            else:
                frame = normalize_scores(source_path, row["name"], source, direction, threshold)
            frame.to_csv(out, index=False)
            outputs.append(out)
        return outputs

    def map(self, name=None, force=False):
        outputs = []
        for row in self.selected(name):
            out = self.config.work_dir / "mapping" / f"{row['name']}.csv"
            if out.exists() and not force:
                outputs.append(out)
                continue
            mappings = build_mapping(self.config.resolve(row["fasta"]), self.config.resolve(row["pdb"]), row["chain"])
            out.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([asdict(item) for item in mappings]).to_csv(out, index=False)
            outputs.append(out)
        return outputs

    def prepare(self, name=None, force=False):
        self.normalize(name=name, force=force)
        self.map(name=name, force=force)
        outputs = []
        mode = self.config.settings.get("online_mode", "mutation_effect")
        if mode == "pairwise_coupling":
            raise ValueError("pairwise coupling tables cannot directly generate single-mutant FoldX jobs")
        saturation = bool(self.config.settings.get("saturation_scan", True))
        for row in self.selected(name):
            normalized = pd.read_csv(self.config.work_dir / "normalized" / f"{row['name']}.csv")
            selected = normalized[normalized["selected"].astype(str).str.lower().isin({"true", "1"})]
            sequence = read_fasta(self.config.resolve(row["fasta"]))
            if saturation:
                candidates = expand_saturation(sorted(set(int(value) for value in selected["sequence_position"])), sequence)
            else:
                candidates = [Candidate(int(item.sequence_position), str(item.wild_type).upper(), str(item.mutant).upper()) for item in selected.itertuples()]
            mapping_frame = pd.read_csv(self.config.work_dir / "mapping" / f"{row['name']}.csv")
            from .models import ResidueMapping
            mapping = {int(item.sequence_position): ResidueMapping(int(item.sequence_position), str(item.sequence_residue), str(item.chain), int(item.structure_position), str(item.structure_residue)) for item in mapping_frame.itertuples()}
            work = self.config.work_dir / "foldx" / row["name"]
            work.mkdir(parents=True, exist_ok=True)
            pdb_target = work / Path(row["pdb"]).name
            if force or not pdb_target.exists():
                shutil.copy2(self.config.resolve(row["pdb"]), pdb_target)
            individual_list = work / "individual_list.txt"
            write_individual_list(candidates, mapping, individual_list)
            pd.DataFrame([asdict(item) for item in candidates]).to_csv(work / "candidates.csv", index=False)
            outputs.append(individual_list)
        return outputs

    def run(self, name=None, force=False):
        self.prepare(name=name, force=force)
        outputs = []
        executable = self.config.resolve(self.config.tools["foldx"])
        for row in self.selected(name):
            work = self.config.work_dir / "foldx" / row["name"]
            runner = FoldXRunner(executable, work, number_of_runs=self.config.settings.get("number_of_runs", 1), out_pdb=self.config.settings.get("out_pdb", False))
            repaired = work / f"{Path(row['pdb']).stem}_Repair.pdb"
            if force or not repaired.exists():
                repaired = runner.repair(Path(row["pdb"]).name)
            dif = runner.build_batch(repaired.name, work / "individual_list.txt")
            outputs.append(dif)
        return outputs

    def report(self, name=None):
        outputs = []
        for row in self.selected(name):
            work = self.config.work_dir / "foldx" / row["name"]
            dif_files = sorted(work.glob("Dif_*.fxout"), key=lambda path: path.stat().st_mtime, reverse=True)
            if not dif_files:
                raise FileNotFoundError(f"no Dif_*.fxout for {row['name']}")
            mutations = [line.strip() for line in (work / "individual_list.txt").read_text(encoding="utf-8").splitlines() if line.strip()]
            frame = parse_dif_fxout(dif_files[0], mutations)
            candidates = pd.read_csv(work / "candidates.csv")
            frame["sequence_position"] = candidates["sequence_position"].values
            frame["wild_type"] = candidates["wild_type"].values
            frame["mutant"] = candidates["mutant"].values
            result_dir = self.config.results_dir / row["name"]
            result_dir.mkdir(parents=True, exist_ok=True)
            frame.to_csv(result_dir / "foldx_results.csv", index=False)
            build_workbook(frame, result_dir / "foldx_results.xlsx")
            build_figures(frame, result_dir)
            outputs.append(result_dir / "foldx_results.csv")
        return outputs

    def all(self, name=None, force=False):
        self.run(name=name, force=force)
        return self.report(name=name)
