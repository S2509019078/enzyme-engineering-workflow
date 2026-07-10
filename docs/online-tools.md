# 在线工具结果导入

导出 EVcouplings/ESM 类结果后放入 `inputs/`，用 `normalize` 转为统一表。至少需要位点列和分数列；可以使用 `position`/`pos`、`wt`、`mutant`/`mutation`、`score`/`ev_score` 等常见列名。分数方向必须在配置中明确，不能根据单个结果自动猜测。

在线网站的版本、队列和服务端参数可能变化，因此仓库保存导出的结果和来源说明，而不声称完全本地重现网站。

