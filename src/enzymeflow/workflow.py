from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import csv
import json
import shutil
import yaml

import pandas as pd

from .cache import build_payload, manifest_valid, write_manifest
from .foldx import FoldXRunner
from .mapping import build_mapping, read_fasta
from .models import Candidate, ResidueMapping
from .mutations import expand_saturation, partition_candidates, write_individual_list
from .normalize import normalize_pairwise_couplings, normalize_scores
from .report import build_figures, build_workbook, parse_dif_fxout


class WorkflowConfig:
    def __init__(self, path: Path):
        self.path = path.resolve()
        if not self.path.exists():
            raise FileNotFoundError(self.path)
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
        self.enzymes = self._load_enzymes(config.enzymes_table) if config.enzymes_table.exists() else []

    @staticmethod
    def _load_enzymes(path: Path):
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
            fields = set(reader.fieldnames or [])
        required = {"name", "fasta", "pdb", "chain", "online_result"}
        if not rows:
            raise ValueError(f"empty enzyme table: {path}")
        missing = required - fields
        if missing:
            raise ValueError(f"missing enzyme columns: {', '.join(sorted(missing))}")
        names = [row["name"] for row in rows]
        if len(names) != len(set(names)):
            raise ValueError("enzyme names must be unique")
        return rows

    def _setting(self, row, key, default):
        value = str(row.get(key, "")).strip()
        return value if value else self.config.settings.get(key, default)

    def _manifest(self, stage: str, name: str) -> Path:
        return self.config.work_dir / "manifests" / stage / f"{name}.json"

    def selected(self, name=None):
        if not self.enzymes:
            raise ValueError(f"no enzyme definitions loaded from {self.config.enzymes_table}")
        rows = self.enzymes if not name else [row for row in self.enzymes if row["name"] == name]
        if name and not rows:
            raise ValueError(f"unknown enzyme: {name}")
        return rows

    def check(self):
        problems = []
        if not self.config.enzymes_table.exists():
            problems.append(f"missing enzyme table: {self.config.enzymes_table}")
            return problems
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
        for row in self.selected(name):
            threshold = float(self._setting(row, "selection_threshold", 0.0))
            direction = self._setting(row, "score_direction", "higher_is_better")
            source = self._setting(row, "online_source", "mutation_effect")
            mode = self._setting(row, "online_mode", "mutation_effect")
            source_path = self.config.resolve(row["online_result"])
            out = self.config.work_dir / "normalized" / f"{row['name']}.csv"
            payload = build_payload({"online_result": source_path}, {"threshold": threshold, "direction": direction, "source": source, "mode": mode})
            manifest = self._manifest("normalize", row["name"])
            if not force and manifest_valid(manifest, payload, [out]):
                outputs.append(out)
                continue
            out.parent.mkdir(parents=True, exist_ok=True)
            frame = normalize_pairwise_couplings(source_path) if mode == "pairwise_coupling" else normalize_scores(source_path, row["name"], source, direction, threshold)
            frame.to_csv(out, index=False)
            write_manifest(manifest, payload, [out])
            outputs.append(out)
        return outputs

    def map(self, name=None, force=False):
        outputs = []
        for row in self.selected(name):
            fasta = self.config.resolve(row["fasta"])
            pdb = self.config.resolve(row["pdb"])
            out = self.config.work_dir / "mapping" / f"{row['name']}.csv"
            payload = build_payload({"fasta": fasta, "pdb": pdb}, {"chain": row["chain"]})
            manifest = self._manifest("mapping", row["name"])
            if not force and manifest_valid(manifest, payload, [out]):
                outputs.append(out)
                continue
            mappings = build_mapping(fasta, pdb, row["chain"])
            out.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([asdict(item) for item in mappings]).to_csv(out, index=False)
            write_manifest(manifest, payload, [out])
            outputs.append(out)
        return outputs

    def prepare(self, name=None, force=False):
        self.normalize(name=name, force=force)
        self.map(name=name, force=force)
        outputs = []
        for row in self.selected(name):
            mode = self._setting(row, "online_mode", "mutation_effect")
            if mode == "pairwise_coupling":
                raise ValueError("pairwise coupling tables cannot directly generate single-mutant FoldX jobs")
            normalized_path = self.config.work_dir / "normalized" / f"{row['name']}.csv"
            mapping_path = self.config.work_dir / "mapping" / f"{row['name']}.csv"
            normalized = pd.read_csv(normalized_path)
            selected = normalized[normalized["selected"].astype(str).str.lower().isin({"true", "1"})]
            if selected.empty:
                raise ValueError(f"no candidates passed the selection threshold for {row['name']}")
            sequence = read_fasta(self.config.resolve(row["fasta"]))
            saturation = str(self._setting(row, "saturation_scan", True)).lower() in {"true", "1", "yes"}
            candidates = expand_saturation(sorted(set(int(value) for value in selected["sequence_position"])), sequence) if saturation else [Candidate(int(item.sequence_position), str(item.wild_type).upper(), str(item.mutant).upper()) for item in selected.itertuples()]
            mapping_frame = pd.read_csv(mapping_path)
            mapping = {int(item.sequence_position): ResidueMapping(int(item.sequence_position), str(item.sequence_residue), str(item.chain), int(item.structure_position), str(item.structure_residue)) for item in mapping_frame.itertuples()}
            valid, rejected = partition_candidates(candidates, mapping)
            work = self.config.work_dir / "foldx" / row["name"]
            work.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([{**asdict(candidate), "reason": reason} for candidate, reason in rejected], columns=["sequence_position", "wild_type", "mutant", "reason"]).to_csv(work / "rejected_candidates.csv", index=False)
            if not valid:
                raise ValueError(f"all candidates were rejected after FASTA-PDB validation for {row['name']}")
            exact_scores = {(int(item.sequence_position), str(item.wild_type).upper(), str(item.mutant).upper()): item for item in selected.itertuples()}
            metadata = []
            for candidate in valid:
                exact = exact_scores.get((candidate.sequence_position, candidate.wild_type, candidate.mutant))
                record = asdict(candidate)
                record.update({"source_score": getattr(exact, "source_score", None) if exact else None, "source": getattr(exact, "source", None) if exact else None, "source_direction": getattr(exact, "source_direction", None) if exact else None, "source_selected_exact_mutation": bool(exact)})
                metadata.append(record)
            batch_size = int(self._setting(row, "batch_size", 200))
            if batch_size < 1:
                raise ValueError("batch_size must be positive")
            expected_batch_names = {f"batch_{index // batch_size + 1:04d}" for index in range(0, len(valid), batch_size)}
            for old in work.glob("batch_*"):
                if old.is_dir() and old.name not in expected_batch_names:
                    shutil.rmtree(old)
            batch_manifest = {"schema_version": 1, "enzyme": row["name"], "batches": []}
            for index in range(0, len(valid), batch_size):
                batch_candidates = valid[index:index + batch_size]
                batch_metadata = metadata[index:index + batch_size]
                batch_dir = work / f"batch_{index // batch_size + 1:04d}"
                batch_dir.mkdir(parents=True, exist_ok=True)
                write_individual_list(batch_candidates, mapping, batch_dir / "individual_list.txt")
                pd.DataFrame(batch_metadata).to_csv(batch_dir / "candidates.csv", index=False)
                batch_manifest["batches"].append(batch_dir.name)
                outputs.append(batch_dir / "individual_list.txt")
            (work / "batch_manifest.json").write_text(json.dumps(batch_manifest, ensure_ascii=False, indent=2), encoding="utf-8")
            source_pdb = self.config.resolve(row["pdb"])
            pdb_target = work / source_pdb.name
            if force or not pdb_target.exists() or source_pdb.read_bytes() != pdb_target.read_bytes():
                shutil.copy2(source_pdb, pdb_target)
        return outputs

    def _batch_names(self, work: Path) -> list[str]:
        manifest = work / "batch_manifest.json"
        if not manifest.exists():
            raise FileNotFoundError(f"missing batch manifest: {manifest}")
        data = json.loads(manifest.read_text(encoding="utf-8"))
        names = data.get("batches", [])
        if not names:
            raise ValueError(f"batch manifest contains no batches: {manifest}")
        return names

    def _batch_result_valid(self, batch_dir: Path, dif: Path, mutations: list[str]) -> bool:
        try:
            frame = parse_dif_fxout(dif, mutations)
        except Exception:
            return False
        return len(frame) == len(mutations)

    def run(self, name=None, force=False):
        self.prepare(name=name, force=force)
        outputs = []
        executable = self.config.resolve(self.config.tools["foldx"])
        for row in self.selected(name):
            work = self.config.work_dir / "foldx" / row["name"]
            number_of_runs = int(self._setting(row, "number_of_runs", 1))
            if number_of_runs != 1:
                raise ValueError("number_of_runs values other than 1 are not yet supported safely")
            root_runner = FoldXRunner(executable, work, number_of_runs=1, out_pdb=False)
            repaired = work / f"{Path(row['pdb']).stem}_Repair.pdb"
            repair_payload = build_payload({"pdb": work / Path(row["pdb"]).name}, {"command": "RepairPDB"}, {"foldx": str(executable)})
            repair_manifest = self._manifest("repair", row["name"])
            if force or not manifest_valid(repair_manifest, repair_payload, [repaired]):
                repaired = root_runner.repair(Path(row["pdb"]).name)
                write_manifest(repair_manifest, repair_payload, [repaired])
            for batch_name in self._batch_names(work):
                batch_dir = work / batch_name
                individual = batch_dir / "individual_list.txt"
                candidates_csv = batch_dir / "candidates.csv"
                mutations = [line.strip() for line in individual.read_text(encoding="utf-8").splitlines() if line.strip()]
                expected = batch_dir / f"Dif_{repaired.stem}.fxout"
                payload = build_payload({"repaired_pdb": repaired, "individual_list": individual, "candidates": candidates_csv}, {"number_of_runs": 1}, {"foldx": str(executable)})
                manifest = self._manifest("build", f"{row['name']}__{batch_name}")
                if not force and manifest_valid(manifest, payload, [expected]) and self._batch_result_valid(batch_dir, expected, mutations):
                    outputs.append(expected)
                    continue
                for stale in batch_dir.glob("Dif_*.fxout"):
                    stale.unlink()
                shutil.copy2(repaired, batch_dir / repaired.name)
                runner = FoldXRunner(executable, batch_dir, number_of_runs=1, out_pdb=False)
                dif = runner.build_batch(repaired.name, individual)
                if not self._batch_result_valid(batch_dir, dif, mutations):
                    raise ValueError(f"incomplete or mismatched FoldX output in {batch_dir}")
                write_manifest(manifest, payload, [dif])
                outputs.append(dif)
        return outputs

    def report(self, name=None):
        outputs = []
        for row in self.selected(name):
            work = self.config.work_dir / "foldx" / row["name"]
            frames = []
            for batch_name in self._batch_names(work):
                batch_dir = work / batch_name
                dif_files = sorted(batch_dir.glob("Dif_*.fxout"))
                if len(dif_files) != 1:
                    raise FileNotFoundError(f"expected one Dif_*.fxout in {batch_dir}, found {len(dif_files)}")
                mutations = [line.strip() for line in (batch_dir / "individual_list.txt").read_text(encoding="utf-8").splitlines() if line.strip()]
                frame = parse_dif_fxout(dif_files[0], mutations)
                candidates = pd.read_csv(batch_dir / "candidates.csv")
                if len(frame) != len(candidates):
                    raise ValueError(f"FoldX row count mismatch in {batch_dir}")
                for column in candidates.columns:
                    frame[column] = candidates[column].values
                frame["batch"] = batch_dir.name
                frames.append(frame)
            combined = pd.concat(frames, ignore_index=True)
            combined["source_rank"] = combined["source_score"].rank(method="min", ascending=False, na_option="bottom") if "source_score" in combined else None
            combined["candidate_priority"] = combined.apply(lambda item: "model_and_stability" if pd.notna(item.get("source_score")) and item["ddg_kcal_mol"] <= -1.0 else ("stability_only" if item["ddg_kcal_mol"] <= -1.0 else "review"), axis=1)
            result_dir = self.config.results_dir / row["name"]
            result_dir.mkdir(parents=True, exist_ok=True)
            combined.to_csv(result_dir / "foldx_results.csv", index=False)
            build_workbook(combined, result_dir / "foldx_results.xlsx")
            build_figures(combined, result_dir)
            outputs.append(result_dir / "foldx_results.csv")
        return outputs

    def all(self, name=None, force=False):
        self.run(name=name, force=force)
        return self.report(name=name)
