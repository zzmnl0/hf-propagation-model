# 真实 IGRF 场下的 O/X 完整反射：结果与误差定位

Task：`AH3D-PH474-S5-REFLECTIGRF-20261103`，2026-11-03。项目 HEAD `4be052c0c69f31a631beb8ed823c13d539036274`。本报告绑定已保存采集与比较，最终接受范围以[实际独立终审](./AH3D-PH474-S5-REFLECTIGRF-20261103-review.md)为准；Worker自检不等于独立验收。

## 本轮完成的科学比较

保留 TRANSNE 的真实 `Ne0(h)*(1+.006*(lon_deg-120))`、球面R=6371km、12MHz/25°/方位0°、原支持域与两端原生常数，只把恒ECEF B替换为冻结IGRF空间场。8条新独立参考、4条AH与4条PH逻辑轨迹均真实返回地面。PH共6次启动/6份raw，无重发：O/X细档各两段原生续接，拼接31,572/33,483节点；O/X粗档分别4,113/971节点。

AH实际使用生产cubic IGRF provider；PH球面B节点对应同一ECEF物理点，经WGS84位置查询同一`ppigrf 2.1.0/IGRF14.shc`、2020-06-01T12:00 UTC，再转ECEF Tesla。90,415网格点及156直接控制共90,571点。两端连续Ne/B插值和原生常数分别保留；这不是两端原生IGRF生成器全域等价，也不是WGS84传播几何已验收。

生产AH、默认、PH包及MATLAB偏好修改均0；原参考和基线轨迹未重跑。工程R0的真实负控制缺口已在同窗口R1补齐并获[实际独立E ACCEPT](./AH3D-PH474-S5-REFLECTIGRF-20261103-engineering-review.md)，[R0失败记录](./AH3D-PH474-S5-REFLECTIGRF-20261103-engineering-review-r0.md)保留。

## 反射参数及整幅磁场替换效应

下表为本端新L2自然端减冻结TRANSNE本端L2自然端。东/北分量使用TX固定ECEF基底，不是经纬度近似换算；变化包括发射点B幅值/方向及空间变化，不能单独称为磁场梯度贡献。

| 端/模式 | 东向落点变化 km | 北向落点变化 km | 群时延变化 μs |
|---|---:|---:|---:|
| AH O | +7.370346 | +36.239828 | +144.599322 |
| PH O | +7.370393 | +36.238190 | +144.592685 |
| AH X | −6.458714 | −26.846251 | −104.498435 |
| PH X | −6.458722 | −26.846764 | −104.500557 |

两端及各自独立参考预测的方向和量级一致；变化量误差仍须对本端参考分开核查，不能称两端变化完全相同。AH变化量对参考的北向误差约+4.72/−15.26mm；PH对应约+0.347/+0.137m。完整几何程/群程/方向变化及其误差见[Worker科学表](./AH3D-PH474-S5-REFLECTIGRF-20261103-worker-report.md)和[原始比较](./AH3D-PH474-S5-REFLECTIGRF-20261103/comparison.json)。本场AH的O−X落点ECEF弦距约139.7429km，群时延差约0.536760ms。

发射及地面均核Ne=gradNe=0后，群传播方向可由波法线推导为真空方向；它不是PH原生群向量遥测。若接收DOA定义指向来波，则取传播向量反号。方向变化误差表是两个无符号夹角的差，不冒称有符号旋转向量。

## 自然端误差与可确认的分解

| L2模式 | AH对本端Q2位置误差 m | PH对本端Q2位置误差 m | AH−PH位置向量模 m |
|---|---:|---:|---:|
| O | 0.0211 | 0.4239 | 0.2089 |
| X | 0.00227 | 0.1963 | 0.6958 |

各自同类真实自然端的逐分量恒等式为 `A−P=(A−RA)−(P−RP)+(RA−RP)`。前两项合为本端数值残差组合；最后一项联合包含原生常数及连续Ne/B定义差异。群程的有符号例子（m）：

| 模式 | 数值组合 | 条件组合 | AH−PH群程 |
|---|---:|---:|---:|
| O | −0.404695 | +0.282673 | −0.122021 |
| X | −0.188530 | +1.038643 | +0.850113 |

完整精度的分解闭合；上表独立舍入可带来末位差。O模式存在数值项与条件项抵消，不能套用前包“条件差主导”的统一结论，也不能用较小AH−PH差直接推断AH自身更准。位置和p向量亦保存有符号分解，范数不能直接相加；每1km的完整共同程量表见[O](./AH3D-PH474-S5-REFLECTIGRF-20261103/O_L2_signed_common_s.npz)、[X](./AH3D-PH474-S5-REFLECTIGRF-20261103/X_L2_signed_common_s.npz)。

**本轮未发现新增AH生产修复依据。** 当前显式SI、FD=0.03125km及L2配置对独立参考的本包自然端、保存节点壳/支持、重建等门通过；不能外推到生产默认、所有节点间状态或一般场。AH两L1自然端严格门通过，但重建/采样门仍FAIL，不迁移为完整通过。

## PH资格和原因边界

PH两条L2自然端对本端参考的有限数值门通过，完整数值及物理资格仍为False：

- O/X原生末状态低于球面约16.03/8.15cm，末点path/state差约同量级，均违反原1e−8km严格门。
- 保存内点最大壳残差约4.44e−8/3.91e−8，超过1e−10；末点壳本身约0/2.22e−16，不能用合格末点掩盖沿程失败。
- PH对本端参考的1km/PCHIP主表最大p差约1.41e−6/1.89e−6，超过1e−6；重建位置影响约0.943/1.032m，超过0.1m。
- AH−PH共同程最大位置/群程虽约0.279/0.105m（O）、0.593/0.231m（X），p差约1.307e−6/1.347e−6仍超原门，不能据位置和群程小于1m宣布整曲线通过。

以上是确认的数值/接口现象。PH投影、积分、输出/事件处理各自的因果份额、D_omega一般符号守卫风险及截止/退化分支仍未确认；不凭这些现象指定单一内部原因，不为追随地下末点修改正确AH地面事件。

## 证据及实际调用链

- AH：`engineering-r1/ah_schedule.py→execute.child→checked_ah.py→run.ah_one`，在真实`model_for`返回处核实际`IGRFGridGeomagneticField3D(cubic)`及`RefractiveIndexAH3DECEF`，由生产`integrate_fixed_path_ah3d_ecef`推进。见[实际调度](./AH3D-PH474-S5-REFLECTIGRF-20261103/engineering-r1/ah-dispatch.jsonl)、[采集代码](./AH3D-PH474-S5-REFLECTIGRF-20261103/run.py)、[生产磁场类](../models/geomagnetic_field_3d.py)、[生产积分器](../models/ray_tracer_ah_3d_ecef.py)。四条记录的Ne/B查询及物理节点支持审计全PASS，越界0，证据在[AH产物](./AH3D-PH474-S5-REFLECTIGRF-20261103/AH/)各case的support-audit。
- PH：实际Root工具调用`launch_ph_r1.ps1→validation.py→launch_ph.ps1→ph_export.m→raytrace_3d_sp`，再由`consume.parse/stitch`核源/12状态/累计量并拼接。见[启动与退出汇总](./AH3D-PH474-S5-REFLECTIGRF-20261103/root-ph-launch/launch-summary.json)及[原始PH](./AH3D-PH474-S5-REFLECTIGRF-20261103/PH/)。六段各抽样最多129保存节点核Ne/B，全PASS，最大B误差约7.58e−20T、Ne相对误差约5.57e−13；这是保存节点抽样，不是全部节点或内部RK查询证明。
- 独立参考：`reference.integrate/dynamics/field/density`保留完整`epsilon_Ne*gradNe+epsilon_B:J_B`，两端连续场分别重构。8参考、四组双档/自然端/壳已被实际E独审复算接受，后续采集未重积分参考。PH球面边缘沿用[既有静态证据](./AH3D-PH474-VERIFY-20261013-mag_field_interp_3d_sp.c:119)，本轮无新REA调用。

逐次PH R1外层stdout未单独落盘，保留为遥测证据缺口；[Root调度见证](./AH3D-PH474-S5-REFLECTIGRF-20261103/root-ph-r1-invocations.json)明确是采集后记录的实际调用命令，不冒称逐次marker回执。[正向live握手](./AH3D-PH474-S5-REFLECTIGRF-20261103/root-r1-live-preflight.log)有真实`R1_PH_DISPATCH_VALIDATED`，请求/raw/runtime/源/接缝与实际6次资源退出回执齐备。没有为补遥测重跑科学。

## 解决清单与总体位置

1. 已有SI常数/实际系数消费/IGRF cubic精度修复继续保留；本包未新增生产修复，不能将对PH的剩余差全部称作AH缺陷。
2. 当前可用成果为限定IGRF真实O/X反射参数及逐量、逐模式的有符号有限归因。PH完整曲线及严格物理资格保持失败；AH粗档重建问题保留，后续按路径/散射用途需要再处理。
3. 未确认项包括一般三维连续场、WGS84真实反射联动、生产默认全域资格、一般D_omega/截止退化及各内部机制独立贡献；它们不是已确认的AH待修bug。
4. [实际Planner建议](./AH3D-PH474-S5-REFLECTIGRF-20261103-next-direction.md)：本包终审接受后，保留当前Ne/IGRF来源并明确同一真实位置映射，下一步仅细化WGS84真实O/X反射，随后S6分层验收与解决清单汇总。不先穷尽PH内部小机制或重跑下降矩阵，本轮不细化/实施下一包。

总体仍在S5的受控组合验证及S6有限归因增量阶段，完整S5与全误差修复尚未完成。正式新SI/旧物理各0/26、历史兼容18/26、F35与旧FAIL/BLOCKED/UNKNOWN分别保留，不迁移计数。最终科学清单、实际终审及closeout绑定另列；活动时钟/独审表不进入冻结科学叶。
