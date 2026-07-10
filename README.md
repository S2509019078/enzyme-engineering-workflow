# EnzymeFlow：可复用酶改造与 FoldX 流程

Windows 优先。新项目只需替换 FASTA、PDB、在线工具导出的 CSV 和配置，不需要修改 Python 代码。

流程把在线序列模型分数与 FoldX ΔΔG 分开保存：前者用于候选优先级，后者用于结构稳定性估计。FoldX 需用户自行安装并在 `config/config.yaml` 指定路径；程序不会下载或分发 FoldX。

## 快速开始

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -e .
Copy-Item config\config.example.yaml config\config.yaml
Copy-Item config\enzymes.example.tsv config\enzymes.tsv
python -m enzymeflow.cli check --config config\config.yaml
```

常用入口：`check` 检查配置，`normalize` 标准化 EV/ESM 导出，`map` 检查 FASTA—PDB 位点对应，`prepare` 生成 19 种突变列表，`run` 执行可断点续跑的 FoldX，`report` 只用已有结果重建表格和热图，`status` 查看状态。

如果状态 CSV 被 Excel 占用，主计算不会终止，会自动写入旁路 JSONL 日志。输出中的 `source_score` 和 `ddg_kcal_mol` 始终是不同列，不能混为同一个分数。

