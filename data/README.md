# 建模数据清单

更新日期：2026-08-11

## 当前模型的数据需求

当前可运行实现位于根目录 `bioconcrete/`，使用方法见 `MODELING.md`。旧的
`perform/bioconcrete_mineralization_model.py` 仅保留为基线原型。

新模型只保留群体级等效动力学、碳酸盐化学、反应传输和材料修复任务，明确排除蛋白序列、菌株
工程、基因回路、突变设计及生物构建流程。`data/processed/model_priors/` 只包含聚合统计。
数据库按下列优先级使用：

1. 公共动力学数据的聚合分布：约束等效消耗速率和半饱和常数的先验范围。
2. 项目材料实验：校准释放、活性、失活和沉淀速率。
3. 碳酸盐与水泥相热力学：计算物种分布、饱和指数和可能生成的固相。
4. 裂缝修复实验：校准闭合率、渗透率和吸水率之间的关系。

## 当前建模结果

最终 28 天先验基线位于 `model_runs/final/`：

| 层级 | 平均裂缝闭合率 | 局部最高闭合率 | 平均导流能力比 | CaCO3 平均质量 | 氨生成 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0D | 2.776% | 2.776% | 0.9190 | 0.0624 kg/m3 | 0 |
| 1D | 1.117% | 2.842% | 0.9672 | 0.1171 kg/m3 | 0 |
| 2D | 0.282% | 3.065% | 0.9917 | 0.0776 kg/m3 | 0 |

这些结果是数据库和文献先验下的未校准预测，不代表最终修复性能，也没有把 80% 目标写入方程。
快速数值验证已通过碳/钙守恒、非负性、无源极限、无氨、时间步和网格收敛检查，报告位于
`model_runs/validation/validation.json`。下一步应使用项目自己的乳酸消耗、CaCO3 质量、裂缝图像、
渗透率或吸水率时间序列进行校准。

## 已下载数据库

| 数据库 | 本地目录 | 已验证内容 | 主要用途 |
| --- | --- | --- | --- |
| BRENDA 2026.1 | `brenda/` | 完整 JSON，676.78 MiB | Km、kcat、pH、温度、底物和生物来源 |
| SABIO-RK | `sabio-rk/raw/` | CA 124 条、脲酶 42 条、LDH 728 条 | 标准化动力学参数、实验条件和速率方程 |
| Rhea | `rhea/raw/tsv/` | EC、UniProt、ChEBI、KEGG 映射和反应 SMILES | 统一反应与化合物标识 |
| UniProt | `uniprot/raw/` | CA 142 条、脲酶 827 条、项目菌蛋白 4688 条 | 蛋白序列、功能、催化反应和辅因子 |
| UniProt 精选集 | `uniprot/selected/` | 8 条项目目标蛋白 | ldh、lutA/B/C 和 4 条碳酸酐酶 |
| BacDive v2 | `bacdive/raw/` | Sutcliffiella cohnii 11 个菌株 | 培养基、温度、pH、需氧性、代谢和安全信息 |
| PHREEQC 3.8.6 | `phreeqc/raw/database/` | 19 个 `.dat` 数据库及水泥相扩展 | 碳酸盐平衡、Calcite、Portlandite、C-S-H/AFm 相 |
| CEMDATA18 | `cemdata/raw/` | 2022-03-31 PHREEQC 数据库 | 完整水泥水化相、Calcite/Aragonite、Portlandite 和 C-S-H |

## 关键发现

- `Bacillus cohnii` 的当前分类名为 `Sutcliffiella cohnii`，UniProt Taxonomy ID 为 33932。
- 项目菌蛋白组包含 `ldh`、`lutA`、`lutB`、`lutC`，支持乳酸利用机制。
- 项目菌包含 4 条带 EC 4.2.1.1 注释的碳酸酐酶候选蛋白。
- BacDive 类型株 1081 是专性需氧、可形成芽孢的菌株；这支持模型保留氧气激活门控。
- BacDive 类型株的脲酶实验结果为阴性，因此脲酶路线不应作为当前项目菌的主矿化机制。
- PHREEQC 的 `Concrete_PHR.dat` 和 `Concrete_PZ.dat` 包含 Portlandite、碳酸盐 AFm、Jennite 和 Tobermorite 类 C-S-H 相。

## 文件使用建议

- 机器学习和统计清洗优先读取 JSON/TSV。
- ODE 速率方程优先读取 SABIO-RK SBML；LDH SBML 已按每页 100 条拆成 8 个有效文件。
- `sabio-rk/quarantine/*.partial` 是失败传输的隔离文件，不得用于分析。
- PHREEQC 通用低离子强度计算先用 `phreeqc.dat`；高离子强度孔隙液考虑 `pitzer.dat`，并按需 `INCLUDE$ Concrete_PZ.dat`。
- 不要在未检查重复相名称时同时加载完整 CEMDATA 和 `Concrete_PHR.dat`/`Concrete_PZ.dat`。

## 仍需实验或文献补齐

标准数据库不能直接提供以下项目特异参数：

- 工程菌在真实混凝土孔隙液中的乳酸钙消耗曲线；
- 微胶囊破裂、释放和扩散折减系数；
- 强碱环境中的菌体失活速率；
- CaCO3 质量、晶型和时间序列；
- 裂缝填充率与渗透率/吸水率之间的标定关系；
- CaCO3/C-S-H 复合填充增强系数。

这些参数应由 wet lab 时间序列和项目材料试验校准，不能由公共数据库中的跨物种数据直接替代。
