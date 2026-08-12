# BioConcrete Multiscale Model

面向自修复混凝土的多尺度物理化学与反应传输模型。

本工程研究一个具体问题：当裂缝进水后，封装修复材料释放，乳酸钙在氧气参与下被消耗，碳酸体系
逐步形成并与钙离子发生矿化时，裂缝内部会在什么位置、以多快的速度产生方解石，以及这些固体沉积
会怎样改变裂缝宽度、孔隙率、渗透率和导流能力。

项目提供从均匀反应体系到真实二维裂缝的三个模型层级、统一参数配置、公共数据库聚合流程、实验
校准接口、敏感性分析、数值验证和可复现的 28 天基线结果。

> **结果性质：**仓库中的数值是公共数据库和文献先验下的未校准模拟，不是湿实验结果，也不是对
> 最终修复性能的承诺。80% 修复率仅是项目评价目标，从未被写入方程或强制拟合。

## 为什么建立这个模型

自修复混凝土同时包含多个相互影响的过程：

1. 水进入裂缝后触发胶囊材料释放。
2. 氧气、底物、pH 和温度共同限制群体级修复活性。
3. 乳酸钙消耗产生钙和无机碳，但部分碳会进入等效生物量。
4. CO2、HCO3- 和 CO3^2- 的比例随 pH 改变。
5. CaCO3 只有在方解石达到过饱和后才会沉淀。
6. 沉淀会阻碍后续扩散，同时减小孔隙率和裂缝导流能力。

只用一个经验“修复率”公式无法解释这些反馈。本模型把反应、传输、矿化和材料性能放入同一计算
框架，使每个参数都有单位、来源和可校准接口，并允许检查质量守恒和网格收敛。

## 工程范围与公开边界

仓库只处理材料与数学建模：

- 胶囊释放及等效群体活性；
- 乳酸钙、氧、钙和总无机碳的反应与传输；
- 碳酸盐平衡、Portlandite 溶解和方解石沉淀；
- 0D、1D 和真实 2D 数值求解；
- 裂缝闭合、孔隙率、渗透率、导流能力和吸水率代理量；
- 参数校准、bootstrap 置信区间、Morris 与 Sobol 敏感性分析。

仓库明确不包含蛋白序列、菌株工程、遗传回路、突变设计、培养方案或生物构建流程。BRENDA 和
SABIO-RK 只被转换为反应类别级的分位数统计，不公开物种、序列或位点信息。

## 模型的整体数据流

```text
公共动力学记录             水泥相热力学/文献参数
       |                            |
       v                            v
聚合参数先验 ----------> 统一配置与参数来源表
                              |
          +-------------------+-------------------+
          |                   |                   |
          v                   v                   v
       0D 反应器          1D 裂缝轴向模型       2D 裂缝模型
          |                   |                   |
          +-------------------+-------------------+
                              |
                              v
          CaCO3/C-S-H 体积、闭合率、孔隙率、渗透率
                              |
                              v
                 实验校准、敏感性与验证报告
```

## 反应模型

### 1. 状态变量

每个空间单元追踪 13 个状态：

| 状态 | 含义 |
| --- | --- |
| `capsule_calcium_lactate_mol_m3` | 胶囊内尚未释放的乳酸钙 |
| `spore_density_rel` | 未激活修复单元的相对密度 |
| `active_density_rel` | 活跃修复单元的相对密度 |
| `lactate_mol_m3` | 已释放乳酸浓度 |
| `oxygen_mol_m3` | 溶解氧浓度 |
| `calcium_mol_m3` | 可反应钙浓度 |
| `inorganic_carbon_mol_m3` | 总无机碳浓度 |
| `hydrated_carbon_mol_m3` | 水合碳酸盐部分 |
| `portlandite_mol_m3` | 有效 Portlandite 储量 |
| `calcite_mol_m3` | 累积方解石沉淀 |
| `csh_volume_fraction` | 胶囊释放的纳米 C-S-H 体积分数 |
| `biomass_carbon_mol_m3` | 为碳守恒设置的等效生物量碳 |
| `ammonium_mol_m3` | 无氨路径诊断量，必须恒为零 |

### 2. 胶囊释放与环境门控

胶囊释放使用受水活度控制的一阶过程：

```text
r_release = k_release * gate_aw * C_capsule
```

修复活性由水活度、氧气、pH 和温度共同决定。阈值使用平滑 Sigmoid 或高斯函数，避免硬开关导致
数值不连续。默认环境为每天湿润 12 小时的间歇喷淋，也支持持续湿润和干燥情景。

### 3. 乳酸与氧限制

等效乳酸消耗率采用双 Monod 方程：

```text
r_L = q_max * B_active
      * L / (K_L + L)
      * O2 / (K_O2 + O2)
      * f_aw * f_pH * f_T
```

模型以守恒反应作为物质计量基础：

```text
Ca(C3H5O3)2 + 6 O2 -> CaCO3 + 5 CO2 + 5 H2O
```

实现时以“每摩尔乳酸”为单位换算，并用 `biomass_carbon_fraction` 扣除进入等效生物量的碳。氧
消耗量随可氧化碳比例同步变化，因此封闭体系可以检查碳和钙守恒。

### 4. 碳酸盐与方解石沉淀

长期模拟默认使用 CO2/HCO3-/CO3^2- 的代数平衡；需要研究早期瞬态时，可将
`carbonate_mode` 改为 `kinetic`，显式计算水合与脱水过程。

方解石饱和度近似为：

```text
Omega_calcite = a_Ca2+ * a_CO3^2- / Ksp
```

只有 `Omega_calcite > 1` 时发生沉淀：

```text
r_calcite = k_calcite * A_s * max(Omega_calcite - 1, 0)^n
```

空间模型使用预生成的地球化学查询表。当前环境没有 PHREEQC 可执行程序，因此查询表采用透明、可
复现的碳酸盐解析近似；元数据会明确记录后端，并检查本地 PHREEQC/CEMDATA 文件是否含 Calcite、
Portlandite、Tobermorite、Jennite 和 C-S-H 等相。它是参数先验，不能替代孔隙液实验。

## 空间模型

### 0D：均匀批量反应器

0D 假定体系空间均匀，用于检查化学计量、质量守恒、反应时间尺度和参数敏感性，也是实验标定的
第一层。局部反应由 SciPy `solve_ivp(method="BDF")` 求解，适合刚性动力学。

### 1D：裂缝口到尖端

1D 模型沿 100 mm 裂缝长度离散，裂缝口连接外界氧和无机碳，尖端采用无通量边界。它回答沉淀能
进入多深、裂缝口附近是否优先封堵，以及胶囊分布如何形成轴向修复差异。

### 2D：真实长度和宽度方向求解

2D 模型直接在裂缝长度和宽度组成的有限体积网格上求解，并在裂缝壁附近放置可复现的离散胶囊
源。它不是把 1D 曲线复制或插值成热图。默认网格为 `15 x 5`，可在配置中提高分辨率。

所有空间层级使用算子分裂：

1. 每个网格单元执行局部 BDF 反应更新；
2. 溶解组分执行隐式有限体积扩散/对流更新；
3. 根据新增固体更新孔隙率、扩散阻力和下一步反应环境。

传输通式为：

```text
d(phi*C_i)/dt = div(phi*D_eff*grad(C_i)) - div(q*C_i) + Q_i
```

## 从矿物体积到修复性能

方解石摩尔数通过摩尔质量和密度转换为体积，再与已释放 C-S-H 体积分数相加：

```text
V_solid = n_calcite * M_calcite / rho_calcite + V_CSH
```

当前几何闭合近似为 `H = V_solid`，并限制在 `[0, 1]`。有效裂缝宽度为
`b = b0 * (1-H)`。材料性能使用两种关系：

- 多孔基体渗透率：相对 Kozeny-Carman 关系；
- 裂缝导流能力：立方律 `T/T0 = (b/b0)^3`。

吸水率代理量取相对渗透率平方根。这里的“代理量”需要用实际吸水实验标定，不能直接当作标准试验
值。

## 28 天基线结果

默认条件为 0.3 mm 宽、100 mm 长裂缝，30 C，pH 11.5，每天湿润 12 小时。结果如下：

| 层级 | 平均闭合率 | 局部最高闭合率 | 平均渗透率比 | 平均导流能力比 | CaCO3 平均质量 | NH4+ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0D | 2.776% | 2.776% | 0.9079 | 0.9190 | 0.0624 kg/m3 | 0 |
| 1D | 1.117% | 2.842% | 0.9626 | 0.9672 | 0.1171 kg/m3 | 0 |
| 2D | 0.282% | 3.065% | 0.9905 | 0.9917 | 0.0776 kg/m3 | 0 |

空间平均值低于局部最高值，是因为胶囊只覆盖部分网格，且氧和底物受扩散限制。更重要的是，这组
结果说明当前先验并不能自动得到目标修复率。需要用项目自己的乳酸消耗、CaCO3 质量、裂缝图像和
渗透/吸水数据校准释放、消耗与沉淀速率。

代表性输出：

### 0D 时间曲线

![0D timecourse](model_runs/final/20260811_231012_0d/timecourse.png)

### 1D 空间剖面

![1D profiles](model_runs/final/20260811_231144_1d/healing_profiles.png)

### 2D 裂缝修复图

![2D healing map](model_runs/final/20260811_231655_2d/healing_map.png)

## 数值验证

快速验证全部通过：

| 检查 | 结果 | 阈值 |
| --- | ---: | ---: |
| 封闭体系碳守恒误差 | `6.06e-16` | `< 0.5%` |
| 封闭体系钙守恒误差 | `1.54e-15` | `< 0.5%` |
| 时间步减半差异 | `1.21e-7` | `< 5%` |
| 1D 粗细网格差异 | `6.40e-5` | `< 5%` |
| 2D 粗细网格差异 | `1.42e-5` | `< 5%` |
| 无修复源时矿化 | `0` | `< 1e-10 mol/m3` |
| 氨生成 | `0` | 必须为 `0` |

快速检查使用缩短的空间模拟筛查离散误差；完整 28 天加密网格验证可通过 `validate --full` 运行，计算
成本明显更高。数值收敛不等于实验有效性，模型仍需要独立湿实验验证。

## 安装

需要 Python 3.8 或更高版本。

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -e .
```

主要依赖为 NumPy、SciPy、pandas 和 matplotlib。

## 快速运行

```powershell
# 查看全部命令
python -m bioconcrete --help

# 生成默认配置
python -m bioconcrete default-config --output model_config.json

# 从本地数据库生成不含序列信息的聚合先验
python -m bioconcrete prepare-data

# 生成碳酸盐地球化学查询表
python -m bioconcrete build-geochem-grid

# 运行三个层级
python -m bioconcrete simulate --level 0d --config model_config.json
python -m bioconcrete simulate --level 1d --config model_config.json
python -m bioconcrete simulate --level 2d --config model_config.json

# 执行数值与物理检查
python -m bioconcrete validate --config model_config.json
```

默认 2D 运行比 0D/1D 显著更慢，因为它会求解完整 28 天湿干循环和每个网格单元的刚性反应。

## 实验校准

校准 CSV 必须包含 `time_d`，并可包含以下任意实测列：

```text
group, crack_width_mm, pH, oxygen_mg_L, lactate_mM,
cfu_mL, caco3_mg, healing_ratio, permeability_ratio, sorptivity
```

运行：

```powershell
python -m bioconcrete calibrate --data experiment.csv --bootstrap 20
```

默认拟合胶囊释放、等效乳酸消耗和方解石沉淀速率，并用 bootstrap 输出参数置信区间。没有实验文件
时，程序不会使用模拟数据进行“伪校准”。

## 敏感性分析

```powershell
python -m bioconcrete sensitivity --config model_config.json --samples 8
```

程序生成 Morris 基本效应和 Sobol 一阶/总效应指数。样本数 8 适合流程测试；正式分析应提高样本
数并报告随机种子、参数范围和收敛情况。

## 仓库结构

```text
bioconcrete/                  核心 Python 包
  config.py                   配置、参数范围与来源类别
  data_pipeline.py            BRENDA/SABIO-RK 聚合统计
  chemistry.py                碳酸盐平衡和地球化学查询表
  model.py                    0D/1D/2D 反应传输求解器
  analysis.py                 校准、bootstrap、Morris、Sobol
  validation.py               守恒、极限情景和收敛检查
  report.py                   CSV、PNG 和 Markdown 报告
  cli.py                      命令行入口
data/processed/               可公开的聚合参数与查询表
model_runs/final/             三个层级的可复现基线结果
model_runs/validation/        数值验证证据
tests/                        自动测试
model_config.json             默认 28 天配置
MODELING.md                   精简技术说明
```

原始数据库、PDF、旧原型、中间运行目录和本地实验文件通过 `.gitignore` 排除。特别是完整 BRENDA
JSON 约 677 MiB，不进入 Git 历史。使用者需要根据各数据库许可自行下载原始数据；仓库只提交聚合
统计和来源说明。

## 当前局限与下一步

1. **缺少项目实验校准。** 当前数值只能作为先验情景，不能用于工程设计或安全决策。
2. **地球化学为解析近似。** 当前未调用 PHREEQC 可执行程序，强离子孔隙液下的活度需要进一步验证。
3. **闭合率关系简化。** 固体体积分数到裂缝几何闭合的映射需要显微图像标定。
4. **C-S-H 仅作为胶囊释放固体。** 模型不会虚构“微生物生成 C-S-H”。
5. **材料力学尚未显式求解。** 当前输出渗透和导流能力，不直接预测抗压/抗折强度恢复。
6. **二维网格为工程折中。** 默认 `15 x 5` 便于复现，正式研究应执行网格加密和全时段验证。

优先实验数据是 0、1、3、7、14、21、28 天的乳酸浓度、CaCO3 质量、裂缝宽度图像、渗透率或
吸水率。完成校准后，应留出独立试验组验证模型，而不是用同一批数据同时拟合和评价。

## 复现性说明

- 默认随机种子为 `2026`，用于生成二维胶囊位置；
- 每个运行目录保存实际配置、状态表、摘要、诊断、图像和报告；
- GitHub Actions 自动执行编译、单元测试和快速物理验证；
- 参数来源区分为模型先验、数据库聚合、文献先验、热力学数据库、项目假设和项目实验。

更紧凑的命令说明见 [MODELING.md](MODELING.md)，数据边界与已有数据说明见
[data/README.md](data/README.md)。
