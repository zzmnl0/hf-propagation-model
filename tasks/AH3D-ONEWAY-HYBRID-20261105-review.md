# AH3D 单程混合模型计划独立审查

| 项目 | 内容 |
|---|---|
| Task ID | AH3D-ONEWAY-HYBRID-20261105 |
| 日期 | 2026-11-05 |
| Reviewer | 实际调用的 /root/scatter_route_review |
| 决定 | ACCEPT，仅覆盖总体计划与 S1 方案 |
| 基线 HEAD | 8fbb5e1210b0a13fb2dcfea34eee86261752bd5e |
| 工作树 | 审查时 tracked 工作树干净 |
| 总体计划 | AH3D-ONEWAY-HYBRID-20261105-plan.md，v1.1 |
| 总体计划 SHA256 | 9CA2196AE3C8F41994DA4B093586672E5ADF5C3A55F6C62DB74B0560C339D8F6 |
| S1 方案 | AH3D-ONEWAY-HYBRID-20261105-S1-plan.md，v1.1 |
| S1 方案 SHA256 | 1B6D748EB80FFBD3B544EB3CDB1962C64ADD1BB2B060DE39F5768670782FFE4D |

本记录由 Root 按独立 Reviewer 实际返回的最终意见保存。ACCEPT 对应规划质量。S1 的执行状态和物理结果由后续实施证据记录。

[总体计划](AH3D-ONEWAY-HYBRID-20261105-plan.md) [S1 详细方案](AH3D-ONEWAY-HYBRID-20261105-S1-plan.md)

## 审查结论

无阻塞发现。

S1 的共同介质快照、功率归一 E/H 端口、双向均匀薄层参考及同基底 S0/ΔS，足以构成进入 S2 受控交接工作的基础。方案已纳入独立圆偏振解析对照、非零时刻介质采样、左右不同外部介质、复相干输入功率检查和背景传播相位。

总体计划明确新混合功能默认关闭并显式启用，在接入阶段分别检查普通 AH 和混合入口。S2 交付平层窗口到 RX 的交接，S3 交付具有横向能量再分布的不规则体散射。

## 检查证据与实施状态

Reviewer 独立核对两份文件哈希、Git 状态及相关源码。现有 berreman_matrix_uniform 调用内部 Maxwell 算子，S1 以无磁 Fresnel–Fabry–Pérot 和平行磁场圆偏振解析式提供独立对照。方案也已处理现有 PE 张量默认系数与 AH 实际系数的差异。

本次审查为只读计划审查。控制例、误差阈值和资源限制属于待实施的验收要求。数值作业、测试与 S1 物理验收尚待执行。

## 延期事项及归属

| 事项 | Reviewer 分类与影响 | 对应阶段 |
|---|---|---|
| AH 任意空间初态接续、有限束和 RX 捕获 | CONFIRMED / MAJOR / DEFERABLE，局部端口与平层算子完成后仍需端到端交接证据 | S2 闭环通过前落实 |
| 两横向维不规则体散射 | CONFIRMED / MAJOR / DEFERABLE，平层算子之外需要横向波矢耦合 | S3 散射结果通过前落实 |
| 发射到接收的绝对功率预算 | CONFIRMED / MAJOR / DEFERABLE，单位端口功率之外需要天线、束输运与接收投影 | S4 绝对功率产品通过前落实 |
| 历史 B1 参考求解器资格 | CONFIRMED / MAJOR / DEFERABLE，原正式细化门仍为资源阻塞 | 引用其正式验收资格前落实，实际采用求解器另有对应阶段证据 |

上述事项均已安排到实际需要它们的阶段，S1 保持当前完整交付范围。

## Planner Feedback

- REQUIRED 保持已冻结的介质、单位、端口规范和三组验收，控制例或阈值发生实质变化时记录原因及验收影响。S2 须证明 AH 接续、接收捕获与窗口替换。
- RECOMMENDED 按当前顺序完成平层交接、横向散射与完整观测量，历史小规模参考结果按其实际范围使用。
- OPTIONAL 本次审查无需增加研究阶段或控制例。
