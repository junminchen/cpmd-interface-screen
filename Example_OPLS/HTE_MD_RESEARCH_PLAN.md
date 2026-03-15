# 高通量实验 + Constant-Potential MD 研究计划

## 1. 研究目标

建立一个闭环研究框架，将高通量实验与 constant-potential 界面 MD 耦合起来，用可解释的界面机制描述符预测并解释锂金属相关电化学表现。

核心目标：

1. 建立 `配方 -> 界面描述符 -> 电化学表现` 的数据链路
2. 验证多电压界面描述符是否优于传统组成特征或 bulk solvation proxy
3. 形成可用于后续配方推荐的监督学习数据集

优先目标变量：

- `coulombic_efficiency_percent`
- `interfacial_resistance_ohm_cm2`
- `capacity_retention_percent`

## 2. 中心假设

与仅使用组成特征或体相溶剂化 proxy 相比，带电压分辨率的界面描述符更接近影响 CE 的真实物理过程，因此对电化学表现应具有更高的解释力和更好的迁移性。

具体假设：

1. 高 CE 配方会表现出更强的界面阴离子参与或更合理的添加剂竞争配位行为
2. 单一 bulk `anion ratio` 不足以表征偏压下真实界面重构
3. 多电压 constant-potential MD 提取的界面代理变量可作为 CE 的高阶 proxy

## 3. 总体工作流

```text
配方库设计
  -> 高通量实验测试
  -> constant-potential MD（多电压）
  -> 界面描述符提取
  -> 数据集汇总
  -> 基线模型 + MD描述符模型 + 融合模型
  -> 误差分析与下一轮配方选择
```

## 4. 配方空间设计

### 4.1 第一阶段范围

建议从 `30-60` 个严格受控配方起步。

控制原则：

- 固定主盐浓度窗口
- 固定主溶剂家族
- 有计划地改变添加剂类别与含量
- 保持测试协议一致

建议变量：

- 盐种：例如 `LiPF6`, `LiFSI`, `LiTFSI`
- 基础溶剂体系：例如 `EC/DMC`, `FEC/EMC`, `DME/DOL`
- 添加剂种类：含氟碳酸酯、醚类、磷酸酯、磺酮类等
- 添加剂含量：低、中、高 3 档

### 4.2 编码规则

每个配方必须有唯一 `formulation_id`。

推荐字段：

- `formulation_id`
- `salt_name`
- `salt_concentration_m`
- `solvent_1`
- `solvent_1_fraction`
- `solvent_2`
- `solvent_2_fraction`
- `additive_1`
- `additive_1_fraction`
- `additive_2`
- `additive_2_fraction`
- `notes`

## 5. 高通量实验设计

### 5.1 必测标签

每个配方至少记录：

- `formulation_id`
- `cell_type`
- `temperature_c`
- `current_density_ma_cm2`
- `areal_capacity_mah_cm2`
- `ce_cycle_1`
- `ce_mean`
- `ce_std`
- `interfacial_resistance_ohm_cm2`
- `capacity_retention_percent`

### 5.2 推荐附加实验信息

建议同步记录：

- 电导率
- 黏度
- Li+ transference number
- 初始过电位
- EIS 等效电路拟合参数
- 代表性样品的 XPS / cryo-TEM / Raman / NMR

### 5.3 协议一致性

高通量实验必须尽可能固定：

- 电极类型
- 隔膜
- 电流密度
- 面容量
- 温度
- 测试截止条件

否则模型会先学到测试条件差异，而不是配方规律。

## 6. MD 设计

### 6.1 第一阶段 MD 目标

对每个配方跑多电压 constant-potential 界面 MD。

建议每个配方至少覆盖：

- `0.5 V`
- `1.0 V`
- `1.5 V`
- `2.0 V`
- 如资源允许再加 `3.0 V`

### 6.2 每个配方的 MD 产物

每个配方在每个电压点输出：

- 拓扑与轨迹
- `electrode_charges.log`
- descriptor JSON

### 6.3 推荐第一阶段描述符

#### 添加剂相关

- `*_interface_enrichment_factor`
- `*_surface_contact_probability`
- `*_first_layer_residence_time_ps`
- `*_li_shell_replacement_fraction_interface`
- `*_li_shell_replacement_fraction_bulk`
- `*_li_shell_replacement_interface_minus_bulk`

#### 阴离子相关

下一步需要新增：

- `interfacial_anion_ratio_at_V`
- `bulk_anion_ratio_at_V`
- `delta_anion_ratio_at_V`
- `anion_enrichment_factor_at_V`
- `li_anion_coordination_number_interface_at_V`
- `li_anion_coordination_number_bulk_at_V`
- `anion_residence_time_near_interfacial_li_at_V`

#### 电响应相关

- `interfacial_capacitance_at_V`
- `d_interfacial_anion_ratio_dV`
- `d_anion_enrichment_dV`
- `d_capacitance_dV`

### 6.4 数据组织方式

建议最终训练表使用宽表结构：

- 一行一个 `formulation_id`
- 列中显式编码电压，例如：
  - `anion_ratio_V0p5`
  - `anion_ratio_V1p0`
  - `anion_ratio_V2p0`
  - `d_anion_ratio_dV_0p5_1p0`

## 7. 建模设计

### 7.1 第一阶段模型对照组

必须同时做三类模型：

1. 组成特征基线
   - 参考 Cui 风格
   - 输入：元素比、分子类别、浓度、含氟量、氧原子比等

2. bulk proxy 模型
   - 输入：导电率、黏度、传统 `anion ratio`、分子静态描述符

3. 界面 MD 描述符模型
   - 输入：本项目的 multi-voltage interfacial descriptors

4. 融合模型
   - 输入：组成特征 + bulk proxy + MD descriptors

### 7.2 推荐算法顺序

先做：

- ridge
- lasso
- elastic net
- random forest
- gradient boosting / xgboost / lightgbm

暂不建议一开始就做：

- attention 模型
- GNN
- raw trajectory neural encoder

原因：

- 第一阶段样本量仍然偏小
- 当前输入是低维、可解释的表格特征
- 先做强基线更利于发表和解释

### 7.3 验证方案

至少做三类验证：

1. 随机划分
2. 按添加剂家族分组划分
3. 按溶剂体系分组划分

评价指标：

- `R2`
- `RMSE`
- `MAE`
- 排序相关系数

## 8. 研究问题与可检验结论

本项目可以回答的核心问题：

1. 多电压界面阴离子 proxy 是否优于 bulk `anion ratio`
2. 添加剂界面富集和停留时间是否能独立提升 CE 预测
3. Li+ 第一壳层中“阴离子 vs 添加剂 vs 溶剂”的竞争关系是否是 CE 的关键驱动项
4. 微分电容或电压响应特征是否能提供额外信息
5. 融合模型是否显著优于 composition-only 模型

## 9. 论文结构建议

### 论文故事线

1. 现有工作主要依赖组成特征或 bulk proxy
2. CE 受偏压界面过程控制，bulk proxy 信息不足
3. 本工作引入多电压 constant-potential MD，构建界面机制描述符
4. 高通量实验提供真实 CE/阻抗/循环标签
5. 界面描述符模型优于传统基线，并揭示关键界面机制

### 图表建议

图 1：

- 项目框架图
- 配方库 + HTE + MD + ML 闭环

图 2：

- 不同电压下界面示意图
- anion/additive/Li+ 局域结构重构

图 3：

- 描述符与 CE 的相关关系
- 单变量机制图

图 4：

- 基线模型 vs MD 模型 vs 融合模型的性能比较

图 5：

- SHAP 或 feature importance
- 指向关键界面机制

图 6：

- 主动推荐出的新配方及实验验证结果

## 10. 下一步执行清单

### 立刻开始

1. 定义实验侧 `formulation_id` 规则
2. 整理高通量实验标签模板 CSV
3. 在描述符脚本中加入显式阴离子 tracking
4. 在描述符脚本中实现单电压 `interfacial_anion_ratio`
5. 设计多电压输出 schema

### 第二阶段

6. 跑首批 `10-20` 个配方的多电压 MD 做 pilot
7. 检查描述符稳定性、方差和共线性
8. 对齐首批实验标签
9. 建立第一版 baseline 模型和 MD 模型

### 第三阶段

10. 扩展到完整配方库
11. 加入融合模型
12. 做主动学习或贝叶斯优化选点

## 11. 当前最合理的近期交付物

如果以最短路径推进，建议先完成这三个交付物：

1. `多电压界面阴离子描述符脚本`
2. `标准化 formulation_id + labels.csv 模板`
3. `首批 pilot 配方的 aggregated descriptor dataset`

这三项完成后，才真正进入监督学习阶段。
