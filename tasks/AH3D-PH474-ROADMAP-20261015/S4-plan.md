# S4 AH 定向修复：物理 SI 默认与吸收系数一致性

## 当前进度 v73

2026-10-29；所属[总体v73](../AH3D-PH474-ROADMAP-20261015-plan.md)。本轮按用户“物理正确性优先、PH为真值参考、继续实施”完成[PHYSALIGN v1](../AH3D-PH474-S4-PHYSALIGN-20261029-plan.md)：普通ECEF AH省略系数采用冻结SI值，ECEF吸收读取射线实际系数；显式覆盖、PH兼容、旧射线回退、共享K_FP和真实Ne保留。生产及工程已获[实际独立ACCEPT](../AH3D-PH474-S4-PHYSALIGN-20261029-engineering-review.md)，三条自然对照已完成，最终接受范围见[独立终审](../AH3D-PH474-S4-PHYSALIGN-20261029-review.md)。

见[总报告](../AH3D-PH474-S4-PHYSALIGN-20261029-report.md)：同真实介质G1/L2中，旧系数→SI默认的AH−PH native终点位置差 **597.093330→0.034623 m**，群程绝对差 **834.370736→0.049698 m**；新默认对独立SI原门PASS。显式PH系数为 **0.057509/0.082251 m**，没有进一步改善，保持物理默认，不为本例拟合参数。三条各一次均退出0、descending ground、完整性和域门PASS，无新PH/参考运行。

**本包两项已证问题已修，整个S4及完整S3仍不宣称完成。** PH原native地面FAIL、严格AH/PH资格和PH_GRID、旧遥测UNKNOWN、独立机制份额、154阻塞及一般域保持；legacy dispersion旧测试配置不兼容列延期。G1碰撞为0，吸收修复的本轮证据是非零碰撞解析/调用回归，不是新增PH吸收轨迹验收。

下一建议转向一般三维B=0受控连续场，随后常磁场O/X，再WGS84/IGRF；新包需固定连续场/梯度和比较语义，先只细化下一项。当前厘米残差按影响插入定向诊断，不先追所有机制份额。新支路不绕过原P3.1/PH_GRID或等同正式P4，尚未细化或实施。

更新前原字节见[历史S4](../AH3D-PH474-S4-PHYSALIGN-20261029/navigation-before-v73/S4-plan.md)，SHA256 `4409A1062672FB371F47247A92FE905313E977BD468ACA765536C86666D9004B`。以下均为历史，旧“默认不变/未授权”不覆盖当前已实施范围，也不恢复旧运行许可。

## 历史方向 v72：建议规划局部 SI 默认迁移

2026-10-29；所属[总体 v72](../AH3D-PH474-ROADMAP-20261015-plan.md)。[S3 有限自然终点比较及误差定位](../AH3D-PH474-S3-SIQUAL-AH2-20261028-comparison-and-actions.md)已交付，最终接受范围见[独立终审](../AH3D-PH474-S3-SIQUAL-AH2-20261028-offline-sc-recovery-final-review.md)。显式 SI 通路、已有 SIACCESS 有限兼容回归及本次 AH−SI 自然终点数值门支持优先处理已证默认常数近似；完整 S3、严格 AH 物理资格及 PH 专属资格尚未成立。

**下一轮只建议细化“AH3D 局部 SI 默认迁移与兼容回归”一个工作包，待用户确认；本轮未细化、批准或实施生产默认修改。** 拟讨论的边界：物理 AH3D 入口省略系数时用 SI；显式覆盖优先，PH 兼容模式及显式旧常数复现保留；审计直接/工厂调用及历史夹具，保留共享 K_FP、真实 Ne 和冻结 P3.1 原契约。独立色散函数的省略参数调用须先查调用者，不盲目扩展整改范围。

不为拟合 PH 离壳末点改 AH 物理，不将当前 B=0 单例推广为一般三维/O-X/WGS84/IGRF 等价；PH FAIL、严格资格 BLOCKED/遥测 UNKNOWN、独立机制份额及 154 局部阻塞保留。后续运行矩阵、预算、代码边界和验收门由下一详细计划明确。

本入口更新前原字节见[历史 S4](../AH3D-PH474-S3-SIQUAL-AH2-20261028/navigation-history-before-v72/S4-plan.md)，SHA256 `BC27D43E5973FE2D6DB1CC6EDA5C1B8476FC91727E845EB75302CC21A4D41DBA`。下列 v26 及更早“最近/下一步/当前”描述仅为历史，不覆盖本方向，不恢复旧窗口或科学额度。

## 历史导航 v26 及前序记录

- 所属总体：[总体主计划v26](../AH3D-PH474-ROADMAP-20261015-plan.md)；状态：`SI_LOCAL_ACCESS_DELIVERY_ACCEPTED / S4_NOT_COMPLETE / PH_ENDPOINT_EVIDENCE_PACKAGE_ACCEPTED / NATEND_R1_E_REPLAN / NATEND_R2_E_REWORK / NATEND_R3_E_REPLAN / NATEND_R4_E0_REWORK / NEXT_EXECUTION_NOT_AUTHORIZED / RETRY_NOT_AUTHORIZED / PRODUCTION_DEFAULT_UNCHANGED`。
- 最近实施：[NATEND-R4来源准入E0最终v1](../AH3D-PH474-S3-NATEND-R4-20261017-plan.md)已获[明确执行授权](../AH3D-PH474-S3-NATEND-R4-20261017-execution-authorization.md)，但唯一尝试为实际独立 **REWORK**；见[执行报告](../AH3D-PH474-S3-NATEND-R4-20261017-execution-report.md)、[实际审查](../AH3D-PH474-S3-NATEND-R4-20261017-execution-review.md)及[外层回执](../AH3D-PH474-S3-NATEND-R4-20261017-execution-receipt.json)。ParseFile通过，初始悬挂子入Job为1/1、Resume后计数门首失；触发计数未落盘，清理复用同门，最终Job活动数及自动成员身份UNKNOWN。已知子exit99、外层工具exit1/no session不等于全Job退出或E0通过。子端双向核验/实际错对象及模拟控制未取得完成证据，显式Python/科学调用0不能认证未知自动成员的类型或数量。

本次已消耗唯一尝试，余窗不重开；完整E、D2及R/S/C仍未放行。当前不新细化返工包。下一方向仅为计数先保存后判断、清理独立记录最终活动数与有界成员身份取证，另行规划和许可；不预设辅助成员豁免或放宽上限，不据工程失败新增AH/PH物理错误。

更新前总体v25及S2/S3/S4原字节见[执行历史](../AH3D-PH474-S3-NATEND-R4-20261017/execution-history/)，[执行导航版本映射](../AH3D-PH474-S3-NATEND-R4-20261017-execution-navigation-version-map.json)仅按原路径及期望SHA将本任务v25引用解至同哈希历史；旧v24及更早映射、冻结数学及原科学FAIL/UNKNOWN保持。历史相对链接依原目录语境解释。
- 更新前总体v24及S2/S3/S4原字节见[R4历史](../AH3D-PH474-S3-NATEND-R4-20261017/history/)，[R4显式版本映射](../AH3D-PH474-S3-NATEND-R4-20261017-navigation-version-map.json)仅将R3执行回执27中的4导航及旧执行映射/收尾审查的4项current-v24引用解至同哈希历史；其余23项原路径、v23及更早解析不变，不改旧证据。历史相对链接依原目录语境解释。
- 前序实施（R3历史时点）为NATEND-R3最终v1，[实际独立REPLAN](../AH3D-PH474-S3-NATEND-R3-20261017-execution-review.md)，见[报告](../AH3D-PH474-S3-NATEND-R3-20261017-execution-report.md)。PID0控制NULL/87通过，来源父唯一0x1000查询NULL/5；具体策略及R2旧原因仍未确认。Job/controller/Python/D2及新增科学均0，R/S/C未开，R3停止时尚未新细化工作包，当时仅建议重规划可核验的启动来源链，不提权、不跳过身份门或使用余时重试。实施前v23四导航见[执行历史](../AH3D-PH474-S3-NATEND-R3-20261017/execution-history/)，[执行版本映射](../AH3D-PH474-S3-NATEND-R3-20261017-execution-navigation-version-map.json)只解planning30的4/26、prelaunch68的4/64及旧映射4项current-v23；其余旧来源/数值不改。
- 前序实施（R2历史时点）：NATEND-R2最终v1已获明确批准，但唯一E包装在父进程OpenProcess返回0后硬停止，见[执行报告](../AH3D-PH474-S3-NATEND-R2-20261017-execution-report.md)、[实际独立REWORK](../AH3D-PH474-S3-NATEND-R2-20261017-execution-review.md)及[外层回执](../AH3D-PH474-S3-NATEND-R2-20261017-execution-receipt.json)。过宽权限请求静态确认，实际因果与最小query成功仍UNKNOWN；Job/controller、Python/guard/D2及新增数学全0，不新增AH/PH物理缺陷结论。下一方向仅最小只读访问修正及LastError留证；重试需另行补充计划和许可，不用余时自动重开，不新增S4修复或默认迁移。
- 更新前总体v21及本入口原字节见[执行历史](../AH3D-PH474-S3-NATEND-R2-20261017/execution-history/)；[执行导航映射](../AH3D-PH474-S3-NATEND-R2-20261017-execution-navigation-version-map.json)仅重定位R2 planning25中的4导航和旧映射4个current绑定，其余21项及全部E/数值来源保持原路径。更早[S4导航v20](../AH3D-PH474-S3-NATEND-R2-20261017/history/S4-plan-before-natend-r2.md)、[R2原映射](../AH3D-PH474-S3-NATEND-R2-20261017-navigation-version-map.json)和[R1映射](../AH3D-PH474-S3-NATEND-R1-20261017-navigation-version-map.json)不改；下文旧记录保留历史语义。
- 最近完成：[AH3D-PH474-S4-SICOEFF-20261017 计划v1](../AH3D-PH474-S4-SICOEFF-20261017-plan.md)已按批准范围实施，见[执行报告](../AH3D-PH474-S4-SICOEFF-20261017-execution-report.md)、[实际独立ACCEPT](../AH3D-PH474-S4-SICOEFF-20261017-execution-review.md)与[最终回执](../AH3D-PH474-S4-SICOEFF-20261017-execution-receipt.json)。
- 最近实际验收：[最小SI局部生产接入与兼容回归](../AH3D-PH474-S4-SIACCESS-20261017-plan.md)及E-R1，见[最终实际报告](../AH3D-PH474-S4-SIACCESS-20261017-execution-r1-report.md)、[实际独立ACCEPT](../AH3D-PH474-S4-SIACCESS-20261017-execution-r1-review.md)与[最终回执](../AH3D-PH474-S4-SIACCESS-20261017-execution-r1-receipt.json)。旧首轮采证失败不追认为成功。
- [SIRESID有限离线归因实际ACCEPT](../AH3D-PH474-S3-SIRESID-20261017-execution-review.md)已区分截面重建与节点组合差；[PHEND R1实际独立ACCEPT](../AH3D-PH474-S3-PHEND-R1-20261017-execution-review.md)又确认了终点输出规则，仍不证明完整自然终点准确性。S3的[NATEND-R2](../AH3D-PH474-S3-NATEND-R2-20261017-plan.md)现为实际E REWORK，D2/R/S/C未运行；旧PLAN ACCEPT不替代实际工程验收。S4只承接独立确证的AH缺陷，不为匹配PH高度置零改物理。本轮仅登记状态，不新增S4补丁、默认迁移或采集。
- 本轮更新前[S4导航原字节副本](../AH3D-PH474-S3-NATEND-20261017/history/S4-plan-before-natend.md) SHA `1F00C9FEDA42C741EE2DAEE2992CAFC702DBA8E10CCF69AD00622E561B6DAB8F`；[新版本映射](../AH3D-PH474-S3-NATEND-20261017-navigation-version-map.json)保存旧PHEND规划绑定。更早PHEND前导航仍见[旧历史](../AH3D-PH474-S3-PHEND-20261017/history/S4-plan-before-phend.md)及[旧版本映射](../AH3D-PH474-S3-PHEND-20261017-navigation-version-map.json)，不更改任何旧执行回执。
- 更新前原字节见[历史S4导航](../AH3D-PH474-S3-SIRESID-20261017/history/S4-plan-before-siresid.md)，SHA `323D3D468A71619915819359AD1A6C02613F13F779FF8463C7577BCEF813FFD8`；旧PLAN审查的3项导航历史绑定见[版本映射](../AH3D-PH474-S3-SIRESID-20261017-navigation-version-map.json)，旧执行回执289项不含这些导航，全部仍按原路径核验。
- 更新前S4原字节见[历史S4入口](../AH3D-PH474-S4-SIACCESS-20261017/history/S4-plan-before-siaccess.md)，SHA `091D75C7903BCDE09915BA8AF689452300EB5D316154642EB6429FED9C799A22`。旧审查仍绑定旧导航，本轮由新PLAN审查另行绑定。
- 原提纲已原字节保存在[历史S4入口](../AH3D-PH474-S4-SICOEFF-20261017/history/S4-plan-before-detail.md)，SHA `61015E4CF2D23E7F8E531350BF6D631364E76CFEE61940785E6619C9C0BE42EF`；原“S2A为当前唯一细化对象”仅保留历史语义。

## 目标

消除有独立物理依据、真实触发和明确 AH 根因的物理或实现错误。针对常数近似的校正应有 SI/CODATA 来源，不以匹配 PH 数值为唯一目标。

## 入口与优先级

[S3A实际独立审查](../AH3D-PH474-S3-B0ATTR-20261017-execution-review.md)已接受限定离线归因及SI候选依据。[SICOEFF实际审查](../AH3D-PH474-S4-SICOEFF-20261017-execution-review.md)随后接受两种B=0介质、两设置、60km前缀的显式SI入口有限样本结果；4/4轨迹及两重建均在原SI参考门内，104记录、61必交及89清单成员完整。S父退出/PID/精确UTC仍UNKNOWN，不认证父exit0；该项列为MINOR/DEFERABLE，下一执行任务应完善采证，不为补历史日志重跑旧采集。

已观察的响应包含系数及相应本端上壳初态变化；未认定全部组合残差来自常数，也未分离积分/投影/事件份额。无需等待无关PH未知全部清零，但不能以拟合PH为修正依据。

SIACCESS已在请求末尾加入默认None的`plasma_coefficient`并由`_ah_model`透传；唯一生产工程控制20次构造通过。科学回归实际完成4条球面B=0、O模、10MHz、各60km工厂SI轨迹，对52唯一冻结参考形成104条两重建记录，全部在原门内；与前包直接SI的同截面同重建位置/群程/动量差最大均0。65科学必交和103清单成员齐备；新工程和科学退出证据不补认SICOEFF或SIACCESS首轮旧UNKNOWN。

这确认最小opt-in接入与有限兼容/科学回归，不改变共享K_FP、config默认或真实Ne，也不把公共fan/P2P/retrace的工程路由控制当作其真实轨迹科学验收。下一包回S3识别残差，不重做已通过的20次构造或4条轨迹；S4其他修复由可识别证据触发。

## 预期产物与细化规则

本子任务已实际交付最小生产调用接入补丁、兼容/无效输入/漏传控制、真实工厂路径有限SI回归、退出与来源证据，并获独立实际接受。统一工厂已支持显式系数，不等于SI已迁移为生产默认，不扩大到全部公共流程或一般连续域验收，也不标记整个S4完成。

本轮只更新状态并细化S3下一包，不创建科学输出目录、执行窗口或生产补丁。新包须另行批准，实际实施且独立审查后才细化当时的下一步；总体六阶段不重排，[S5](./S5-plan.md)与[S6](./S6-plan.md)继续仅保留提纲。S3未识别机制、PH两个自然末点FAIL、完整自然终点NOT_ESTABLISHED、一般三维/O-X、WGS84/IGRF、154阻塞及正式P3.1 0/26全部保留。

