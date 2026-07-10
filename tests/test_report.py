from pathlib import Path
from enzymeflow.report import classify_ddg, parse_average_fxout

def test_parse_foldx_average_and_classify(tmp_path):
    p=tmp_path/"Average_x.fxout"
    p.write_text("Pdb\tTotal Energy\nAA1V\t-2.5\nAA1G\t1.2\n", encoding="utf-8")
    df=parse_average_fxout(p)
    assert list(df["mutation"]) == ["AA1V","AA1G"]
    assert classify_ddg(-2.5,-1.0,1.0)=="stabilizing"

