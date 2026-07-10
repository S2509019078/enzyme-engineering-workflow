from pathlib import Path
import subprocess
from .state import StatusJournal

class FoldXRunner:
    def __init__(self, executable, root, journal=None):
        self.executable=str(executable); self.root=Path(root); self.journal=journal or StatusJournal(self.root/"status.jsonl")
    def _run(self, args, key):
        out=self.root/"logs"/f"{key}.log"; out.parent.mkdir(parents=True, exist_ok=True)
        p=subprocess.run([self.executable,*map(str,args)], cwd=self.root, text=True, capture_output=True, check=False)
        out.write_text(p.stdout+"\n"+p.stderr, encoding="utf-8", errors="replace")
        self.journal.safe_status({"task":key,"status":"success" if p.returncode==0 else "failed","returncode":p.returncode})
        return p.returncode
    def repair(self, pdb_name): return self._run(["--command=RepairPDB","--pdb="+str(pdb_name)], "repair_"+Path(pdb_name).stem)
    def build_batch(self, pdb_name, individual_list): return self._run(["--command=BuildModel","--pdb="+str(pdb_name),"--mutant-file="+str(individual_list)], "build_"+Path(individual_list).stem)

