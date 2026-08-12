# BioConcrete Multiscale Model

面向自修复混凝土的多尺度反应传输模型。当前 `v0.2.0` 聚焦于材料建模、公开数据证据链与可复现的不确定性分析。生物侧只通过匿名的群体级等效活性参数进入方程，不保存、推断或推荐生物设计信息。

> **证据边界：**团队目前没有自己的湿实验时间序列。仓库不会把公开数据称为团队实验，也不会用合成数据替代科学校准。由于当前运行环境访问 Zenodo 返回 HTTP 403，公开校准和外部验证结果尚未生成；代码、数据清单、追溯规则与手动下载路径已经实现。

## 工程在做什么

模型将以下过程连接为一个可检验的计算链：

1. 修复材料随湿润条件释放；
2. 匿名等效活性受水活度、氧、pH、温度、变化率和持续时间控制；
3. 乳酸钙、氧、钙和无机碳按守恒关系变化；
4. 动态 pH 由总碱度和 DIC 电荷平衡求解；
5. 方解石仅在过饱和时沉淀，Portlandite 作为碱度和钙的材料储库；
6. 0D、1D 和真实 2D 有限体积模型计算反应与裂缝传输；
7. 固体体积通过壁面沉积厚度映射到裂缝闭合、渗透率和导流能力。

主输出不再把固体体积分数直接当作修复率，而是区分：

```text
solid_fill_fraction
wall_deposition_thickness_mm
crack_closure_ratio
```

其中 `healing_ratio` 只作为旧文件的弃用别名保留。

## 证据等级

| 等级 | 用途 | 当前状态 |
| --- | --- | --- |
| `project experiment` | 团队自己的实验 | 空，0 条 |
| `public calibration data` | 参数估计和内部按试件留出测试 | 数据集 A 等待手动下载 |
| `external validation data` | 冻结参数后的独立验证 | 数据集 B 等待手动下载 |
| `measurement error data` | 裂缝宽度测量噪声 | krkCMd 访问受限 |
| `literature prior` | 参数范围，不作为项目实测 | 已支持 |
| `model prediction` | 待验证的设计建议 | 已预注册，尚未正式批处理 |

数据清单见 [data/public/DATASETS.yml](data/public/DATASETS.yml)，手动下载路径见 [data/public/MANUAL_DOWNLOADS.md](data/public/MANUAL_DOWNLOADS.md)。每条规范化记录保留 DOI 对应的数据集编号、原始文件、工作表和行号。

## 模型层级

- `0d`：均匀批量体系，用于守恒、校准和情景分析。
- `1d`：裂缝口至尖端的轴向反应传输，尖端无通量。
- `2d`：长度和宽度方向均直接求解，并使用未解析厚度计算总质量。

局部刚性反应使用 `scipy.integrate.solve_ivp(method="BDF")`；传输使用隐式稀疏有限体积更新。默认支持每天 12 小时湿润，也支持持续湿润与干燥。

## 动态环境接口

公开参数仅表示跨尺度等效活性：

```text
effective_kcat_s, effective_km_mol_m3, active_unit_concentration,
activity_multiplier, response_delay_h, basal_leak_fraction,
activation_duration_h
```

模型同时输出环境信号、激活状态、激活记忆、累计活性、真激活指数与误激活指数。这些量不能解释为某个具体生物构件的实测性能。

## 安装

Python 3.8 及以上：

```powershell
python -m pip install -r requirements-lock.txt
python -m pip install -e . --no-build-isolation
```

基础模拟不要求 SALib；`formal-sensitivity` 需要锁定的 `SALib==1.4.7`。

## 常用命令

```powershell
python -m bioconcrete default-config --output model_config.json
python -m bioconcrete prepare-data
python -m bioconcrete audit-units

python -m bioconcrete fetch-public-data --manifest data/public/DATASETS.yml
python -m bioconcrete prepare-public-data --dataset transet_18clsu02
python -m bioconcrete prepare-public-data --dataset marine_external

python -m bioconcrete build-phreeqc-grid
python -m bioconcrete compare-geochem-backends

python -m bioconcrete simulate --level 0d --config model_config.json
python -m bioconcrete simulate --level 1d --config model_config.json
python -m bioconcrete simulate --level 2d --config model_config.json

python -m bioconcrete calibrate-public --train data/public/derived/transet_18clsu02/observations.csv --bootstrap 100 --profile-points 9
python -m bioconcrete validate-external --dataset data/public/derived/marine_external/observations.csv --frozen-run model_runs/public_calibration
python -m bioconcrete fit-measurement-error --dataset PATH_TO_KRKCMD_TABLE.csv

python -m bioconcrete formal-sensitivity --samples 1024
python -m bioconcrete design-matrix --preregister PREREGISTERED_SCENARIOS.yml
python -m bioconcrete evidence-report --run model_runs/public_calibration
python -m bioconcrete validate
```

旧 `sensitivity` 命令仅为历史复现入口，不属于正式 Sobol 证据。正式分析使用 SALib，输出 `S1`、`ST`、bootstrap 95% 置信区间和未收敛标记，不会裁剪越界指标。

## 校准与盲测保护

- 数据集 A 按 `specimen_id` 分组划分，禁止同一试件的时间点跨训练/测试集合。
- 每个自由参数至少需要两个独立观测约束。
- 拟合配置保存 SHA-256；外部验证发现配置被修改时立即终止。
- 外部验证代码路径不含优化器，并同时报告零矿化与一阶经验基线。
- krkCMd 只拟合裂缝宽度测量误差，不参与反应动力学拟合。

## 地球化学声明

当前本机只有 PHREEQC 源码包和数据库，没有可执行后端。因此生成表明确标记为 `analytical_surrogate`，`compare-geochem-backends` 返回 `claim_allowed: false`。在真正由 PHREEQC 生成网格并完成饱和指数误差和相分类一致率比较前，项目不得使用“PHREEQC 耦合”表述。

## 预注册设计矩阵

[PREREGISTERED_SCENARIOS.yml](PREREGISTERED_SCENARIOS.yml) 在外部验证前固定了 1,728 个组合：裂缝宽度、剂量、湿润时间、等效活性、响应延迟和基础泄漏。输出采用 Pareto 排序并标记为：

```text
recommended
robust_alternative
not_recommended
```

所有结果必须写作“公开数据支持的待验证预测”，不能描述成团队已经采用或证实的方案。

## 已验证的软件性质

当前 10 项单元测试全部通过。快速物理验收结果：

- 封闭体系碳守恒相对误差：`4.55e-16`；
- 封闭体系钙守恒相对误差：`2.66e-15`；
- 无修复源时矿化量：`0`；
- 氨诊断量：恒为 `0`；
- 时间步、1D 网格和 2D 网格快速收敛差异：均小于 `5%`。

这些结果证明数值实现通过当前测试，不证明材料模型已经获得实验有效性。完整说明见 [WIKI_MODEL.md](WIKI_MODEL.md)，模型公式与实现边界见 [MODELING.md](MODELING.md)。

## 许可证与引用

代码使用 MIT License。引用信息见 [CITATION.cff](CITATION.cff)。
