import argparse
from pathlib import Path

from .wizard import interactive_wizard
from .workflow import EnzymeWorkflow, WorkflowConfig


def main(argv=None):
    parser = argparse.ArgumentParser(prog="enzymeflow", description="可复用酶改造与 FoldX 工作流")
    parser.add_argument("command", choices=["wizard", "check", "normalize", "map", "prepare", "run", "report", "status", "all"])
    parser.add_argument("--config", type=Path, default=Path("config/config.yaml"))
    parser.add_argument("--runs-dir", type=Path, default=Path("runs"), help="wizard 创建独立运行目录的位置")
    parser.add_argument("--enzyme")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    if args.command == "wizard":
        interactive_wizard(base_dir=args.runs_dir)
        return 0

    config = WorkflowConfig(args.config)
    workflow = EnzymeWorkflow(config)

    if args.command == "check":
        problems = workflow.check()
        if problems:
            for problem in problems:
                print(problem)
            return 2
        print("environment OK")
        return 0
    if args.command == "status":
        status_files = sorted(config.work_dir.glob("foldx/*/status.jsonl"))
        if not status_files:
            print("no tasks")
            return 0
        for path in status_files:
            print(f"[{path.parent.name}]")
            print(path.read_text(encoding="utf-8"))
        return 0
    if args.command == "normalize":
        workflow.normalize(name=args.enzyme, force=args.force)
    elif args.command == "map":
        workflow.map(name=args.enzyme, force=args.force)
    elif args.command == "prepare":
        workflow.prepare(name=args.enzyme, force=args.force)
    elif args.command == "run":
        workflow.run(name=args.enzyme, force=args.force)
    elif args.command == "report":
        workflow.report(name=args.enzyme)
    elif args.command == "all":
        workflow.all(name=args.enzyme, force=args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
