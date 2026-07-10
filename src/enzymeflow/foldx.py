from pathlib import Path
import subprocess

from .state import StatusJournal


class FoldXRunner:
    def __init__(self, executable, root, journal=None, number_of_runs=1, out_pdb=False):
        self.executable = str(executable)
        self.root = Path(root)
        self.number_of_runs = int(number_of_runs)
        self.out_pdb = bool(out_pdb)
        self.journal = journal or StatusJournal(self.root / "status.jsonl")

    def _run(self, args, key):
        log = self.root / "logs" / f"{key}.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        process = subprocess.run([self.executable, *map(str, args)], cwd=self.root, text=True, capture_output=True, check=False)
        log.write_text(process.stdout + ("\n" if process.stdout and process.stderr else "") + process.stderr, encoding="utf-8", errors="replace")
        self.journal.safe_status({"task": key, "status": "success" if process.returncode == 0 else "failed", "returncode": process.returncode, "log": str(log)})
        if process.returncode != 0:
            raise RuntimeError(f"FoldX failed for {key}; see {log}")
        return log

    def repair(self, pdb_name):
        pdb_name = Path(pdb_name).name
        self._run(["--command=RepairPDB", f"--pdb={pdb_name}"], "repair_" + Path(pdb_name).stem)
        repaired = self.root / f"{Path(pdb_name).stem}_Repair.pdb"
        if not repaired.exists():
            raise FileNotFoundError(f"FoldX repaired PDB not found: {repaired}")
        return repaired

    def build_batch(self, pdb_name, individual_list):
        pdb_name = Path(pdb_name).name
        mutant_file = Path(individual_list).name
        args = ["--command=BuildModel", f"--pdb={pdb_name}", f"--mutant-file={mutant_file}", f"--numberOfRuns={self.number_of_runs}", f"--out-pdb={1 if self.out_pdb else 0}"]
        self._run(args, "build_" + Path(mutant_file).stem)
        dif = self.root / f"Dif_{Path(pdb_name).stem}.fxout"
        if not dif.exists():
            candidates = sorted(self.root.glob("Dif_*.fxout"), key=lambda path: path.stat().st_mtime, reverse=True)
            if not candidates:
                raise FileNotFoundError("FoldX Dif_*.fxout output not found")
            dif = candidates[0]
        return dif
