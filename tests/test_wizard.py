from datetime import datetime

import yaml

from enzymeflow.wizard import create_run_directory, detect_pdb_chains, write_project


def test_detect_pdb_chains_preserves_order(tmp_path):
    pdb = tmp_path / "test.pdb"
    pdb.write_text(
        "ATOM      1  N   ALA B   1       0.000   0.000   0.000  1.00 20.00           N\n"
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      3  C   ALA B   1       2.000   0.000   0.000  1.00 20.00           C\n",
        encoding="utf-8",
    )
    assert detect_pdb_chains(pdb) == ["B", "A"]


def test_create_run_directory_is_isolated(tmp_path):
    now = datetime(2026, 7, 11, 9, 30, 0)
    first = create_run_directory(tmp_path, "enzyme project", now=now)
    second = create_run_directory(tmp_path, "enzyme project", now=now)
    assert first != second
    assert (first / "inputs" / "fasta").is_dir()
    assert (first / "inputs" / "pdb").is_dir()
    assert (first / "inputs" / "online").is_dir()
    assert second.name.endswith("_01")


def test_write_project_creates_relative_paths_and_settings(tmp_path):
    run_dir = create_run_directory(tmp_path, "demo", now=datetime(2026, 7, 11, 9, 30, 0))
    row = {
        "name": "demo",
        "fasta": "inputs/fasta/demo.fasta",
        "pdb": "inputs/pdb/demo.pdb",
        "chain": "A",
        "online_result": "inputs/online/demo.csv",
        "score_direction": "lower_is_better",
        "selection_threshold": -1.0,
        "saturation_scan": "false",
    }
    config_path = write_project(run_dir, row, foldx_path="foldx")
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    assert config["paths"]["work"] == "work"
    assert config["tools"]["foldx"] == "foldx"
    assert (run_dir / "config" / "enzymes.tsv").exists()
    assert (run_dir / "RUN_INFO.txt").exists()
