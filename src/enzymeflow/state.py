from pathlib import Path
import json, os

class StatusJournal:
    def __init__(self, path): self.path=Path(path); self.backup=self.path.with_suffix(".jsonl")
    def append(self, record):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        try:
            with self.path.open("a", encoding="utf-8") as fh: fh.write(json.dumps(record, ensure_ascii=False)+"\n")
        except OSError:
            with self.backup.open("a", encoding="utf-8") as fh: fh.write(json.dumps(record, ensure_ascii=False)+"\n")
    def safe_status(self, record): self.append(record)

