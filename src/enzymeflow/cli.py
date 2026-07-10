import argparse
from pathlib import Path

def main(argv=None):
    ap=argparse.ArgumentParser(prog="enzymeflow", description="可复用酶改造与 FoldX 流程")
    ap.add_argument("command", choices=["check","normalize","map","prepare","run","report","status","all"])
    ap.add_argument("--config", type=Path, default=Path("config/config.yaml"))
    ap.add_argument("--enzyme")
    ap.add_argument("--force", action="store_true")
    args=ap.parse_args(argv)
    if args.command == "check":
        print(f"config: {args.config} {'OK' if args.config.exists() else 'MISSING'}")
        print("FoldX is supplied by the user and is never bundled.")
        return 0 if args.config.exists() else 2
    if args.command == "status":
        p=Path("work/status.csv"); print(p.read_text(encoding="utf-8") if p.exists() else "no tasks")
        return 0
    print(f"stage={args.command} enzyme={args.enzyme or 'all'} force={args.force}")
    return 0

if __name__ == "__main__": raise SystemExit(main())

