<div align="center">
  <img src="docs/assets/bioconcrete-model-icon.png" alt="BioConcrete Multiscale Model 标志" width="150">
  <h1>BioConcrete Multiscale Model</h1>
  <p><strong>连接修复剂激活、矿物沉积、裂缝闭合与实验设计的多尺度反应传输模型。</strong></p>
  <p>
    <a href="https://github.com/Zwhua/bioconcrete-multiscale-model/releases/tag/v0.5.1"><img src="https://img.shields.io/badge/version-v0.5.1-176B87" alt="v0.5.1 版本"></a>
    <a href="https://github.com/Zwhua/bioconcrete-multiscale-model/actions/workflows/tests.yml"><img src="https://github.com/Zwhua/bioconcrete-multiscale-model/actions/workflows/tests.yml/badge.svg" alt="模型测试"></a>
    <a href="https://github.com/Zwhua/bioconcrete-multiscale-model/actions/workflows/tests.yml"><img src="https://img.shields.io/badge/tests-82%20passing-24745E" alt="82 项测试通过"></a>
    <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/Python-3.8%2B-3776AB" alt="Python 3.8 或更高版本"></a>
    <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-E69F00" alt="MIT 许可证"></a>
  </p>
  <p><a href="README.md">English</a> | 简体中文</p>
  <p>
    <a href="#关键结果"><strong>关键结果</strong></a> ·
    <a href="#证据状态"><strong>证据状态</strong></a> ·
    <a href="#快速开始"><strong>快速开始</strong></a> ·
    <a href="MODELING.md"><strong>技术模型</strong></a>
  </p>
</div>

> [!IMPORTANT]
> **证据状态：** 未校准的机理模型；团队湿实验时间序列数据为 **0 行**。下列数值均为模型输出，不是实验观测。

<table align="center">
  <tr>
    <td align="center"><strong>2.0819%</strong><br><sub>28 天裂缝闭合</sub></td>
    <td align="center"><strong>99.92%</strong><br><sub>闭合来自 C-S-H 先验</sub></td>
    <td align="center"><strong>0 行</strong><br><sub>团队湿实验时间序列</sub></td>
  </tr>
</table>

## 项目概览

| | 摘要 |
| --- | --- |
| **问题** | 预测自修复混凝土裂缝内释放、反应、传输、CaCO3 沉淀和 C-S-H 填充如何共同影响修复。 |
| **模型** | 耦合 0D 动力学、1D 有限体积传输和真实 2D 反应传输，并显式核算材料库存与裂缝几何。 |
| **关键发现** | 默认未校准基线在 28 天达到 **2.0819%** 闭合，其中约 **99.92%** 来自设定的 C-S-H 载荷。 |
| **下一项实验** | 优先比较完整体系、无 C-S-H 条件和非生物 C-S-H 对照，而不是先追求更高的总体活性。 |

这些是**前瞻性模型结论**，尚未被团队实验或公开数据校准证实。

> [!TIP]
> **模型指导的决策：** 优先比较完整体系、无 C-S-H 条件和非生物 C-S-H 对照。由于当前基线主要受库存和几何影响，而不是活性影响，提高总体活性的实验优先级较低。

## 模型价值

模型通过群体尺度的聚合参数连接生物设计和材料性能，但不声称解析某个具体构建体。它区分 CaCO3 沉积与非生物 C-S-H 填充，跟踪碳、钙和完整修复剂库存，并把壁面沉积体积转换为具有明确几何定义的裂缝闭合率。

项目实现了敏感性、实际可识别性、反事实瓶颈和数值 D-optimal 实验设计接口，用于判断下一步最值得测量的变量。模型保留不支持原假设的结果：在当前先验下，只提高总体活性很难显著改善闭合，除非同时改变库存、传输、沉积或几何约束。

## 从模型到决策

<p align="center">
  <img src="model_runs/v0.5.0/figures/figure01_model_to_decision.png" alt="从模型到决策的结构图" width="100%">
</p>

<p align="center"><sub><em>模型结构示意图，不是实验结果。匿名设计类别被映射为可测参数、多尺度材料响应和可证伪实验。</em></sub></p>

## 关键结果

| 结果 | 当前模型输出 | 解释 |
| --- | ---: | --- |
| 28 天裂缝闭合 | **2.0819%** | 未校准基线 |
| C-S-H 贡献 | **2.08023 个百分点** | 当前先验中的主导贡献 |
| 方解石贡献 | **0.00168 个百分点** | 当前先验下贡献较小 |
| 相对导流能力 | **0.9388** | 立方律模型输出 |
| 团队湿实验数据 | **0 行** | 尚未实验验证 |

> [!NOTE]
> **主要结论：** 当前先验配置并非主要受活性限制。若库存、传输、沉积或几何约束不变，仅提高总体活性不太可能带来显著闭合收益。约 2.08% 的结果不是性能承诺，80% 目标也没有被硬编码。

## 证据状态

| 证据组成 | 状态 | 允许的表述 |
| --- | --- | --- |
| 数值守恒 | 已完成 | 通过数值验证 |
| 物理极限与可复现性测试 | 已完成 | 通过数值验证 |
| 未校准 0D 基线 | 已完成 | 仅称为模型输出 |
| 1024 基础样本 Sobol 分析 | 待完成 | 仅有接口，不声明结果 |
| 完整 1,728 情景矩阵 | 待完成 | 已预注册，不声明结果 |
| 公开数据校准 | 等待数据 | 不声明已校准 |
| 独立外部评估 | 等待数据 | 不声明已验证 |
| 团队湿实验时间序列 | 0 行 | 不作实验结论 |

仓库中的 V5 release manifest 是**初始化记录**，不是已完成的正式运行：其中记录 commit `4266ef9`、`git_worktree_dirty: true` 和 `status: initialized`。从 clean commit 生成正式 manifest 仍是待办；仓库不会篡改该来源记录来制造完成状态。

## 快速开始

```bash
git clone https://github.com/Zwhua/bioconcrete-multiscale-model.git
cd bioconcrete-multiscale-model
python -m pip install -r requirements-lock.txt
python -m pip install -e . --no-build-isolation

python -m bioconcrete default-config --output config-v0.5.1.json
python -m bioconcrete simulate --level 0d --config config-v0.5.1.json
python -m bioconcrete validate --config config-v0.5.1.json
```

方程与高级流程见 [MODELING.md](MODELING.md)，公开数据整理见 [data/public/README.md](data/public/README.md)，面向 iGEM 的模型叙事见 [WIKI_MODEL.md](WIKI_MODEL.md)。

## 已实现能力

| 能力 | 实现 | 输出 |
| --- | --- | --- |
| 0D 动力学 | BDF 反应求解器 | 时间过程与守恒量 |
| 1D 传输 | 隐式有限体积模型 | 轴向分布 |
| 2D 传输 | 真实长度-宽度求解器 | 空间分布图 |
| 化学 | 动态电荷平衡；解析地球化学代理 | pH 与碳酸盐状态 |
| 不确定性 | 可续算的先验预测流程 | 先验区间 |
| 实验设计 | 反事实与数值 D-optimal 接口 | 前瞻性决策表 |
| 可复现性 | 随机种子、哈希、冻结配置和 manifest | 可审计运行 |

当前快速验证通过碳钙守恒、非负性、无来源对照、时间步收敛及 1D/2D 网格收敛。本地完整测试集为 **82 项且全部通过**。数值验证不能代替实验有效性。

## 生物设计连接

| 设计类别 | 模型参数 | 所需测量 |
| --- | --- | --- |
| 表面定位 | 有效表面保留/活性单元浓度 | 表面关联的保留活性 |
| 连接区可及性 | 活性倍数 | 匹配条件下的相对活性 |
| 候选酶活性 | 等效动力学先验 | 聚合动力学与稳定性 |
| 微胶囊 | 释放速率与基础泄漏 | 释放曲线和破裂前损失 |
| C-S-H 载荷 | 初始非生物填充 | 非生物 C-S-H 对照 |

这些只是参数映射，不代表模型已测得具体生物实现的真实性能。该接口不保存序列、突变、载体或构建流程。

## 仓库结构

```text
bioconcrete/       核心模型与分析包
data/              先验、数据结构和公开数据清单
model_runs/        按版本保存的分析产物
tests/             数值与证据边界测试
MODELING.md        技术模型文档
WIKI_MODEL.md      面向 iGEM 的叙事文档
```

## 局限与路线图

1. **公开校准数据尚未整理完成。** 动力学和材料响应参数不能称为已校准。
2. **独立外部评估尚未完成。** 当前没有实验 MAE、RMSE、R2、AIC/AICc 或覆盖率结果。
3. **团队湿实验时间序列不可用。** 设计建议仍是可证伪计划，不是已完成的 DBTL 闭环。
4. **PHREEQC 交叉检查待完成。** 当前使用 `analytical_surrogate`，不声称已实现 PHREEQC 耦合。
5. **正式长时分析待完成。** 1024 基础样本 Sobol、完整 1,728 情景矩阵和第 2-8 张结果图有意保持缺失。
6. **裂缝几何经过简化。** 模型未解析完整三维粗糙度、结构力学或现场尺度安全性。

## 引用与许可

当前版本：[v0.5.1](https://github.com/Zwhua/bioconcrete-multiscale-model/releases/tag/v0.5.1)。请按照 [CITATION.cff](CITATION.cff) 引用软件；项目未声明论文或 DOI。代码采用 [MIT License](LICENSE) 发布。
