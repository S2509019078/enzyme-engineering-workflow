# EnzymeFlow：可复用酶改造与 FoldX 工作流

EnzymeFlow 用于把在线序列模型或共进化分析结果，转换为可追踪的候选突变、FASTA—PDB 位点映射、FoldX 批量任务和结构稳定性报告。新项目只需替换 FASTA、PDB、在线工具导出的 CSV 和配置，不需要修改核心 Python 代码。

## 推荐：交互式新建项目

运行：

```powershell
enzymeflow wizard
```

向导会依次询问项目名称、FASTA、PDB、在线预测结果、用于映射的 PDB 链、分数方向、筛选阈值、是否做19种饱和突变以及 FoldX 路径。程序会自动复制输入文件、检测 PDB 链并生成独立运行目录：

```text
runs/20260711_093000_project_name/
  config/
  inputs/fasta/
  inputs/pdb/
  inputs/online/
  work/
  results/
  logs/
  RUN_INFO.txt
```

不同运行的输入、中间文件、FoldX结果和日志不会混用。若 FoldX 已加入系统 PATH，向导输入 `foldx` 即可；若提供明确文件路径，则会保存为绝对路径。

## 工作流边界

在线模型分数和 FoldX ΔΔG 是两类不同证据：前者用于候选位点或候选突变优先级，后者用于估计突变对蛋白结构稳定性的影响。程序始终分别保存 `source_score` 和 `ddg_kcal_mol`，不会把二者混为同一个指标。

支持两类在线结果：

- `mutation_effect`：单突变效应表，至少包含位点或 `A123V` 形式的突变，以及分数；可直接筛选后生成 FoldX 单突变任务。
- `pairwise_coupling`：`i,j,cn` 一类的成对耦联表，只用于共进化位点对分析，不能直接转换成单突变 FoldX 任务。

## 安装

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -e ".[test]"
```

也可以跳过向导，复制示例配置：

```powershell
Copy-Item config\config.example.yaml config\config.yaml
Copy-Item config\enzymes.example.tsv config\enzymes.tsv
```

FoldX 由用户自行安装，并在 `config/config.yaml` 的 `tools.foldx` 中填写路径。仓库不包含 FoldX、本地许可证、密码、Token 或个人路径。

## 输入

`config/enzymes.tsv` 每行对应一个酶，至少包含：

- `name`：项目名；
- `fasta`：目标蛋白 FASTA；
- `pdb`：本地 PDB；
- `chain`：需要映射的结构链；
- `online_result`：在线模型导出的 CSV。

FASTA 与 PDB 不要求编号一致。程序使用序列比对建立真实位点映射，可处理信号肽、前肽、缺失残基和结构编号偏移等常见情况。进入 FoldX 前会检查候选野生型残基、FASTA 残基和 PDB 残基是否一致。

## 运行

```powershell
enzymeflow check --config config\config.yaml
enzymeflow all --config config\config.yaml
```

也可分阶段运行：

```powershell
enzymeflow normalize --config config\config.yaml
enzymeflow map --config config\config.yaml
enzymeflow prepare --config config\config.yaml
enzymeflow run --config config\config.yaml
enzymeflow report --config config\config.yaml
enzymeflow status --config config\config.yaml
```

只处理某个酶：

```powershell
enzymeflow all --config config\config.yaml --enzyme example
```

`--force` 用于覆盖已有中间结果。

## FoldX 任务

程序生成标准 `individual_list.txt`，格式为：

```text
QA213K;
```

即：野生型氨基酸 + 链 + 结构残基号 + 突变氨基酸 + 分号。

默认执行：

- `RepairPDB`
- `BuildModel`
- `--numberOfRuns=1`
- `--out-pdb=0`

最终 ΔΔG 从 `Dif_*.fxout` 的 `Total Energy` 列读取，不从 `Average_*.fxout` 读取。

## 输出

- `work/normalized/`：标准化在线结果；
- `work/mapping/`：FASTA—PDB 位点映射；
- `work/foldx/<enzyme>/batch_*/individual_list.txt`：FoldX 批次突变列表；
- `work/foldx/<enzyme>/status.jsonl`：运行状态；
- `results/<enzyme>/foldx_results.csv`：逐突变 ΔΔG 与在线模型分数；
- `results/<enzyme>/foldx_results.xlsx`：可筛选工作表；
- `results/<enzyme>/foldx_ddg_heatmap.png`：突变氨基酸 × 序列位点二维热图。

## 测试

```powershell
python -m compileall src tests
pytest -q
```

GitHub Actions 会在 Python 3.10、3.11 和 3.12 上运行。单元测试覆盖交互式项目生成、在线结果标准化、FASTA—PDB 比对、FoldX 突变语法、`Dif_*.fxout` 解析和报告逻辑。真实 FoldX 执行仍需本机安装 FoldX。
