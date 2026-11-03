# S5 当前：限定反射结果保留，生产/科研统一L2已ACCEPT（v101）

2026-11-03，所属[总体v101](../AH3D-PH474-ROADMAP-20261015-plan.md)。用户直接要求后续均采用L2，[迁移结果](../AH3D-PH474-L2-DEFAULT-20261103-report.md)及[独立终审](../AH3D-PH474-L2-DEFAULT-20261103-review.md)已ACCEPT；[候选53叶](../AH3D-PH474-L2-DEFAULT-20261103/candidate-manifest.json)SHA `11209981EBAEE5EDD549F46D10DFB763C142B0463D944EB917E4AF6DF2817BB0`，终审SHA `061194417D350DC6A52FCE6DB36B75E792BD1A101E8E7CD452833520A55F3C4B`，[收尾](../AH3D-PH474-L2-DEFAULT-20261103/closeout.json)绑定当前导航。

后续正式AH生产/科研的模型、fan、P2P及实际最终回算均L2，旧档只底层诊断/历史，独立参考精度不降。atol_p=1e−12按k0转物理k；初始步长至多.002km、最大.5km，自适应可更小，名义FD=.03125km且细结构可缩小。P2P同档复用本次L2轨迹，未增加重复积分。既有场景/物理门和独立参考不因默认迁移改变。

新8条短轨+physalign2条、路由/CLI及标签补改定点检查通过；完整8轨之后pipeline显示fallback及测试断言发生两项字节变化，定点验证有界，未假称完整8轨使用最终两文件字节。原科学数据和清单保留，活动源升级通过before快照映射解释；见[版本核验](../AH3D-PH474-L2-DEFAULT-20261103/historical-version-resolution.json)。未重跑旧WGS自然轨迹、PH、正式矩阵或完整P2P，不迁移旧通过数。

S5的既有WGS84/IGRF限定结论继续按原版本解释；PH严格FAIL、机制UNKNOWN及正式0/26等状态保持。后续方向为S6集中增量汇总，不自动扩新控制场。[v100原字节](../AH3D-PH474-L2-DEFAULT-20261103/baseline/S5-v100-before-L2-default.md)已保存，下方为历史。

# 历史 v100：WGS84真实地面O/X反射限定交付已ACCEPT

2026-11-03，所属[总体v100](../AH3D-PH474-ROADMAP-20261015-plan.md)。[WGS84结果](../AH3D-PH474-S5-REFLECTWGS-20261103-report.md)、[实际独立终审](../AH3D-PH474-S5-REFLECTWGS-20261103-review.md)及[收尾](../AH3D-PH474-S5-REFLECTWGS-20261103/closeout.json)已接受限定交付；科学528叶SHA `463DDB039720077D6F7B6713E4C1C30767B227B3DB4C302DC640186ACA89B37F`，审查SHA `8128788FB82AD95ECB69EF1F4072E4DE60B0FAE6523B474BDC86784F6D261012`。8参考/4AH/6PH启动构成4PH逻辑轨迹，均自然GROUND，IGRF新520点，生产/默认/PH包修改0。

AH L2本端自然位置9.183/12.726mm、数值/重建门通过，L1重建FAIL；独审全部9451保存节点壳通过但细档余量很小，未保存RK仍不约束。AH−PH自然端位置5.185649/.695356m；群程6.292159/.944471m由数值组合−.431556/−.229509m与联合条件+6.723715/+1.173980m闭合，不拆机制份额。

PH条件参考不是物理参考通过，严格地面/path-state、壳和共同程/重建FAIL，完整资格False；模板支持另行通过。快逆吻合仅已核非终点样本，终点仍有.155173/.119794m差异，内部RK及机制份额UNKNOWN。新−旧包含真实地面、起点、当地水平与场采样联合变化，不单归曲率或宣称全面WGS84/IGRF等价。正式各0/26及历史18/26保持。

高精度主路线与必要粗档保留，旧V2/历史脚本按[逻辑归档](../AH3D-PH474-AH-CLEANUP-20261103-report.md)退出活跃入口。下一方向为S6集中汇总，不自动追加科学控制或改生产物理；本包后停止。[v99原字节](../AH3D-PH474-S5-REFLECTWGS-20261103/baseline/S5-v99-before-science-closeout.md)已保存，下方为历史。

# 历史 v99：工程及参考进行中

2026-11-03，所属[总体v99](../AH3D-PH474-ROADMAP-20261015-plan.md)。[REFLECTWGS最终v1](../AH3D-PH474-S5-REFLECTWGS-20261103-plan.md)获[实际PLAN ACCEPT](../AH3D-PH474-S5-REFLECTWGS-20261103-plan-review.md)，用户选择真实WGS84站及低层真空延续后按[授权](../AH3D-PH474-S5-REFLECTWGS-20261103-execution-authorization.md)实施。新工程/参考准备中，不授新WGS科学通过；实际工程/参考独审后条件采4AH及最多16PH启动/4逻辑轨迹。显式SI与高精度L2保留，粗档仅必要诊断，核心/默认/PH包修改0。

同ECEF原球面Ne和真实IGRF、椭球地面/发射、原球面顶与实际插值模板支持分别解释。PH完整3D裁剪、快逆/native梯度须实际工程门；条件参考和物理壳/接口严格资格不混同。旧球面对照变化含地面/TX/方向联合因素，不单称曲率。保留所有旧FAIL/UNKNOWN和正式0/26，不先追机制百分比。

[AH脚本整理已限定ACCEPT](../AH3D-PH474-AH-CLEANUP-20261103-review.md)，旧V2原位逻辑归档，冻结路径不移动，不阻塞本包。实际WGS终审后仅建议S6汇总方向，本轮不细化下一包。[v98原字节](../AH3D-PH474-S5-REFLECTWGS-20261103/baseline/S5-v98.md)保存，下方为历史。

# S5 历史：冻结IGRF真实O/X反射限定ACCEPT（v98）

2026-11-03，所属[总体v98](../AH3D-PH474-ROADMAP-20261015-plan.md)。[REFLECTIGRF科学报告](../AH3D-PH474-S5-REFLECTIGRF-20261103-report.md)、[实际独立终审](../AH3D-PH474-S5-REFLECTIGRF-20261103-review.md)及[收尾](../AH3D-PH474-S5-REFLECTIGRF-20261103/closeout.json)已交付。8参考、4AH、6PH启动/6raw形成4逻辑轨迹，均真实GROUND；本包PH O_L1已取得，旧TRANSNE缺口不回填。IGRF共90,571点，实际生产cubic与PH同物理点球面采样分别定义连续场；传播几何仍为球面。生产/默认/PH包/MATLAB偏好改动0。

L2自然端AH对本端参考位置误差O/X约21.1/2.27mm，PH约0.4239/0.1963m；AH−PH约0.2089/0.6958m。AH当前显式FD、保存节点严格/支持及重建门通过，L1重建FAIL。PH地面/path-state、壳、曲线p和重建门仍FAIL，完整数值/物理资格False；自然端有限比较不能替代整曲线验收。O群程数值项与条件项抵消，逐模式有符号分解保留，不能统一称条件主导或拆为独立机制百分比。

本包无新增AH修复依据。下一步仅按[方向反馈](../AH3D-PH474-S5-REFLECTIGRF-20261103-next-direction.md)滚动细化WGS84真实O/X反射，再做[S6分层汇总](./S6-plan.md)；不重复下降矩阵或为逐次stdout缺口补跑。整场IGRF替换含初始B变化，不等于独立梯度贡献或原生IGRF全域等价。PH内部机制、一般截止/退化与全域资格未知继续保留，正式新SI/旧物理各0/26、F35及旧缺口保持。完整S5未完成，本包后停止，下一包未细化/实施。[v97原字节](../AH3D-PH474-S5-REFLECTIGRF-20261103/baseline/S5-v97-before-closeout.md)保存，下方为历史。

# S5 历史：冻结IGRF真实O/X反射工程与参考进行中（v97）

2026-11-03，所属[总体v97](../AH3D-PH474-ROADMAP-20261015-plan.md)。[REFLECTIGRF最终v1](../AH3D-PH474-S5-REFLECTIGRF-20261103-plan.md)已[实际PLAN ACCEPT](../AH3D-PH474-S5-REFLECTIGRF-20261103-plan-review.md)。保留当前Ne/球面/F25，只换真实IGRF空间B；先核两端真实坐标/连续B与完整梯度及8参考，实际E接受后条件采4AH/最多16PH启动，尚无新科学通过。

生产/默认/PH包修改计划0，原阈值/旧严格失败及未知保持；整体场替换效应不解释成独立梯度贡献。下一WGS84包本轮不细化或实施。[v96原字节](../AH3D-PH474-S5-REFLECTIGRF-20261103/baseline/S5-v96.md)保存，下方为历史。

## 历史：单一横向Ne真实O/X反射限定ACCEPT（v96）

2026-11-02，所属[总体v96](../AH3D-PH474-ROADMAP-20261015-plan.md)。[TRANSNE报告](../AH3D-PH474-S5-TRANSNE-20261102-report.md)及[实际独立终审](../AH3D-PH474-S5-TRANSNE-20261102-review.md)已交付。8参考、4AH和9PH启动尝试（8raw/1 O_L1启动失败）及一次离线分析已完成，生产/默认/PH包改动0，旧数据未重跑。4AH及O_L2/X_L1/X_L2真GROUND，两条PH细档自然端齐；O粗档未采集、PH O双档收敛仍未确认。

主L2新−旧基线TX东向位移：AH O/X −4.375988/−3.691147km，PH −4.376008/−3.691145km；时延变化AH+.304803/+.042881μs、PH+.304841/+.042947μs。AH L2自然端位置对本端参考16.228/17.807mm，完整查询支持与已保存节点的相应门通过；AH L1重建FAIL、L2壳有限余量及显式FD配置限定保持。

AH−PH落点距离1.777291/.184132m，组合数值残差.016935/.005066m，联合条件差1.793959/.188284m，逐分量闭合，不代表各机制独立份额，无新增AH生产修复依据。PH重建.340/.362m、壳与地面/末点严格门FAIL，完整资格False，不能因限定ACCEPT迁移为通过。完整数字见[比较](../AH3D-PH474-S5-TRANSNE-20261102/comparison.json)及[收尾](../AH3D-PH474-S5-TRANSNE-20261102/closeout.json)；495科学叶/362工程叶/旧267及15保护均绑定。

下一方向仅为[冻结IGRF空间变B下真实反射，再WGS84反射](../AH3D-PH474-S5-TRANSNE-20261102-next-direction.md)，保留当前Ne和其余条件、复用旧JOINT桥接，不重复下降控制。O_L1缺口延期；完整一般S5、独立机制、F35与正式0/26保持。本包终审后停止，下一包未细化/实施，结果交[S6](./S6-plan.md)有限汇总。[v95原字节](../AH3D-PH474-S5-TRANSNE-20261102/baseline/S5-v95-before-closeout.md)已保存，下方为历史。

# S5 历史：真实反射链路单一横向Ne因素，工程与参考进行中（v95）

2026-11-02，所属[总体v95](../AH3D-PH474-ROADMAP-20261015-plan.md)。[TRANSNE最终v1](../AH3D-PH474-S5-TRANSNE-20261102-plan.md)已[实际独立PLAN ACCEPT](../AH3D-PH474-S5-TRANSNE-20261102-plan-review.md)，依据[本轮继续实施授权](../AH3D-PH474-S5-TRANSNE-20261102-execution-authorization.md)进入工程/参考，尚无新场科学通过。

当前仅在F25真实反射锚点增加正的线性经度因子，保持恒定ECEF B、原发射/域/正确SI及PH原生常数。横向梯度、实际PH连续场/裁剪和8新参考先验证，实际E独审后条件采4AH及最多16PH启动含续接，120分钟共享执行、独审另计。自然端有效性独立于整曲线重建门，不降低原阈值或覆盖PH严格FAIL。

复用已测JOINT桥接，不重复下降控制；本包生产/PH包修改0，未知及正式0/26保持。旧证据/15保护按开始身份核验，实际终审后停止；[v94原字节](../AH3D-PH474-S5-TRANSNE-20261102/baseline/S5-v94.md)保存，下方为历史。

# S5 历史：磁化反射锚点PH自然端已补齐，受限比较完成（v94）

2026-11-02，所属[总体v94](../AH3D-PH474-ROADMAP-20261015-plan.md)。[PH原生续接报告](../AH3D-PH474-S5-REFLECTOX-RESTART-20261102-report.md)及[实际独立终审](../AH3D-PH474-S5-REFLECTOX-RESTART-20261102-review.md)已交付：新增PH六次、O/X各三次，六个接缝12项状态逐值相同，两条细档均真实GROUND。新增AH/参考求解及生产/PH包改动均0；[收尾](../AH3D-PH474-S5-REFLECTOX-RESTART-20261102/closeout.json)绑定267项科学叶、30项工程叶、原164项及15保护文件。

自然端AH−PH位置距离1.773851/0.191752m，组合数值残差0.019972/0.008122m、组合条件差1.793337/0.184450m；有符号分解闭合，主要表现为原生常数与连续Ne插值的组合条件差，不是独立机制比例。本轮没有新增AH修复依据。前包AH细档O/X的有限物理资格保持。

PH完整精度/物理验收仍未通过：重建位置敏感性0.340/0.362m超过原0.1m，壳与地面半径/末点path-state严格门FAIL，`limited_numerical_qualification=False`。本轮ACCEPT仅覆盖原生续接、自然端数据及如实报告的受限比较；AH粗档采样失败、F35及正式新SI/旧物理各0/26均保持。

下一步仅滚动细化[单一Ne横向因素＋恒定ECEF B方向](../AH3D-PH474-S5-REFLECTOX-RESTART-20261102-next-direction.md)，以当前地面发射、电离层转向及真实返回作锚，验证三维偏折和反射参数。复用[JOINT坐标/IGRF桥接](../AH3D-PH474-S5-JOINT-20261101-report.md)，不重复其既定初态下降GROUND控制。连续场及重建能力是必要前置，不能放宽门以解释较大效应。一般三维反射、WGS84/IGRF全域及独立机制仍未确认；本包终审后停止，下一包未细化/实施。

当前实测交[S6](./S6-plan.md)有限汇总；六阶段不重构。[更新前v85原字节](../AH3D-PH474-S5-REFLECTOX-RESTART-20261102/baseline/S5-v85-before-reflection-closeout.md)保存，下方为历史，不覆盖当前v94。

# S5 历史：联合传播有限验收完成，已测覆盖转S6（v85）

2026-11-01：[JOINT报告及解决清单](../AH3D-PH474-S5-JOINT-20261101-report.md)、[实际独立终审](../AH3D-PH474-S5-JOINT-20261101-review.md)、[累计问题台账](../AH3D-PH474-S5-JOINT-20261101-cumulative-issues.md)。非均匀Ne＋冻结IGRF的6新入口、12参考、4AH／4成功PH均已完成，146项manifest和15旧dirty无漂移，本包无生产变更。PH启动实际5次（首失败1，恢复四目标成功）；首次插件内部原因UNKNOWN。

AH O/X两档物理门和真实下降ground通过，PH本条件数值通过，但严格地面／path-state／壳四FAIL且O/X均含大量非末失壳；H只有单档，PH/H条件参考壳全FAIL。L2/PCHIP共同位置O20.795／X16.312mm，群程4.923／6.993mm；主要为连续采样／坐标／原生J组合条件差，不是各机制独立份额。

按[下一方向](../AH3D-PH474-S5-JOINT-20261101-next-direction.md)，停止自动扩增S5控制场，已测有限范围交[S6](./S6-plan.md)集中归因、修复版本验收和正式路线建议。本轮不细化或执行S6，完整S5一般域未验收；PH严格资格／内部RK、一般Domega／截止退化及全部旧FAIL/BLOCKED/UNKNOWN、P3.1 **0/26**保持。[v84原字节](../AH3D-PH474-S5-JOINT-20261101/navigation-before-v85/)已保存，下方为历史。

# S5 历史实施：非均匀Ne＋IGRF联合传播（v84）

2026-11-01：[JOINT最终v1](../AH3D-PH474-S5-JOINT-20261101-plan.md)已获[PLAN ACCEPT](../AH3D-PH474-S5-JOINT-20261101-plan-review.md)，工程开始，本包尚无新轨迹。复用IGRF冻结数组／已修Provider／旧controls与18局部；新Ne全链及6入口、12参考和激发实际独立ACCEPT后条件采4AH4PH。15旧dirty及旧证据保护，新IGRF生成0。

本包接受后停止自动增域，转S6有限汇总方向，不称一般三维全面等价，不在本轮实施S6。PH严格／X内点壳FAIL及全部旧阻塞保持；当前实际资格以后续终审为准。[v83原字节](../AH3D-PH474-S5-JOINT-20261101/navigation-before/)已存，下方均为历史。
# S5 当前：IGRF 有限 ACCEPT；联合传播后转 S6（v83）

2026-11-01：[总报告与解决清单](../AH3D-PH474-S5-IGRF-20261101-report.md)、[实际独立终审](../AH3D-PH474-S5-IGRF-20261101-review.md)。本包已完成18联合局部／12参考／4AH4PH／6生成器点；生产cubic精度缺陷已最小修复并回归，14旧dirty不变，126项manifest独立无漂移。AH有限物理通过，PH条件数值通过，PH严格地面／path-state／壳FAIL保留，X两档有内点壳FAIL。完整S5未完成。

L2/PCHIP共同截面位置O20.934／X18.253mm、群程O5.191／X8.387mm；主要为连续采样／坐标／原生J组合条件差，有符号闭合0，不给独立积分／投影／事件比例。AH自然端位置对参考3.30–6.45微米。旧solver轨迹收益未测，不把节点修复收益冒称轨迹成对收益。

按[下一方向](../AH3D-PH474-S5-IGRF-20261101-next-direction.md)，下一次仅细化同冻结IGRF＋已验非均匀Ne完整O/X组合；独立验收后停止自动扩大控制矩阵，转S6作有限范围归因／修复汇总。下一包尚未细化或执行，完整一般域、PH严格资格、PH_GRID／G1／154／UNKNOWN／P3.1 0/26保持。

v82原导航见[历史](../AH3D-PH474-S5-IGRF-20261101/navigation-before-v83/)，下方为v82及更早历史。
# S5 当前实施：IGRF 桥接与有限 O/X（v82）

2026-11-01：[IGRF最终v1](../AH3D-PH474-S5-IGRF-20261101-plan.md)已获[PLAN ACCEPT](../AH3D-PH474-S5-IGRF-20261101-plan-review.md)，开始工程；当前无本包新轨迹。实际场采样／导数、Ne适配／联合局部及独立参考接受后条件进入4AH／4PH，详见计划及其必要门。

最小cubic精度修复须先重现／定位和独立核具体diff；不以任务场替代错误生产采样。共享节点不等于连续场相同；PH条件重放和严格物理资格分别报告。上一VARB有限ACCEPT及全部旧FAIL／BLOCKED／UNKNOWN保持，完整S5未完成。上一导航原字节见[历史](../AH3D-PH474-S5-IGRF-20261101/navigation-before/)，下方为v81及更早历史。
# S5 分层全程验证：受控平滑变 B 有限验收完成

所属[总体v81](../AH3D-PH474-ROADMAP-20261015-plan.md)。2026-10-30状态`VARB_LIMITED_DIAGNOSTIC_ACCEPT / AH_OX_PHYSICS_PASS / PH_CONDITIONAL_NUMERICAL_PASS / PH_STRICT_NATIVE_FAIL / PH_X_ENDPOINT_SHELL_FAIL`。[VARB总报告与解决清单](../AH3D-PH474-S5-VARB-20261030-report.md)、[实际独立终审](../AH3D-PH474-S5-VARB-20261030-review.md)、[Worker报告](../AH3D-PH474-S5-VARB-20261030-worker-report.md)。18局部、10参考、4AH／4PH及终审完成，科学重跑／参考重算0，112项manifest独立匹配，14份既有dirty不变；完整S5未完成。

单WGS84、真实均匀Ne、Bz经度线性2μT/deg的散度0场，新增磁场空间项及实际Ne/B／经度Jacobian均通过本包必要门。AH自然端位置对本端参考约0.99–2.24μm，L2/PCHIP原生AH−PH共同截面位置O0.295224mm／X0.409592mm、群程O0.895107mm／X0.021344mm；原生自然端位置差O2.352557mm／X40.178126mm。

PH严格地面／path-state四条FAIL，X两档方向壳额外报告门FAIL仅发生末点，非末点≤1.074e−12。原门、原报告及前WGS科学REPLAN不改；有限ACCEPT不含PH-X末点壳或严格接口通过，不作事件／积分／投影份额归因。本包未新增AH生产修复，不调正确AH去拟合PH。

下一轮只细化[真实IGRF生成／采样及生产适配桥接](../AH3D-PH474-S5-VARB-20261030-next-direction.md)，先实际连续场／导数及联合非均匀Ne局部语义，再条件有限O/X轨迹；不先另开完整组合矩阵或PH末点旁支。下一包尚未细化／实施。总体六阶段、PH_GRID、G1／154／UNKNOWN、一般截止退化和P3.1 0/26保持。v80原字节见[导航历史](../AH3D-PH474-S5-VARB-20261030/navigation-before-v81/)，下方为历史，不覆盖v81。

## v80 WGS84 已完成（历史）

所属[总体v80](../AH3D-PH474-ROADMAP-20261015-plan.md)。2026-10-30状态 `WGS84_COLLECTION_COMPLETE / SCIENTIFIC_REPLAN / AH_LIMITED_PHYSICS_PASS / PH_NUMERICAL_PASS / PH_STRICT_NATIVE_FAIL`。[WGS84总报告与解决清单](../AH3D-PH474-S5-WGS84-20261030-report.md)、[独立终审及坐标增量复核](../AH3D-PH474-S5-WGS84-20261030-review.md)、[Worker报告](../AH3D-PH474-S5-WGS84-20261030-worker-report.md)。4AH＋4PH采集、比较及终审已完成，不重跑；完整S5未完成。

同一ECEF位置、波法线及均匀真实Ne/B，B0＋常B O，独立vg直线及各自球／椭球交点参考，新增参考ODE0。AH四条严格物理／地面PASS；PH四条本端数值PASS，但原1e−8km地面及path/state严格门FAIL，总判REPLAN。16组两重建条件组合通过，105项manifest独立核验无漂移。WGS−球面约3.3km自然端差为不同真实地面的解析预期。

PH首点−25.46072mm报告高度偏差已由旧REA `wgs84_xyz2llh_est_` 公式重放并独立复核；非末点最大重放高度差≤3.7251e−12km。原生末点仍在地面下，积分／事件独立份额未确认；不能把正确AH地面事件调向PH地下末点，也不能断言该逆变换只影响输出。本包无新增AH生产修复依据。

下一轮只细化[受控平滑空间变B全程物理方向](../AH3D-PH474-S5-WGS84-20261030-next-direction.md)，同包先核实际Ne/B在共同ECEF位置的场值、方向及相关导数映射；之后再进入真实IGRF／统一工厂。无需先另开整包非均匀WGS或穷尽PH毫米事件份额。本包任务级Provider与PH显式state不代表统一WGS工厂、IGRF、默认角度入口或非均匀完整传播已通过。下一包尚未细化或实施；总体六阶段及原失败／阻塞保持。下方v79及更早状态为历史，不覆盖v80。

## v79 TOPSUP 已验收（历史）

所属[总体v79](../AH3D-PH474-ROADMAP-20261015-plan.md)。2026-10-29状态 `TOPSUP_LIMITED_DELIVERY_ACCEPT / O_FULL_NATURAL_PHYSICS_PASS / PH_STRICT_NATIVE_FAIL`。[TOPSUP 总报告与解决清单](../AH3D-PH474-S5-TOPSUP-20261029-report.md)、[独立终审](../AH3D-PH474-S5-TOPSUP-20261029-review.md)、[下一方向](../AH3D-PH474-S5-TOPSUP-20261029-next-direction.md)。

同解析三维T1、40μT恒ECEF B只扩top460→600/inner450→590，原公式460尺度及五数组旧51层bitwise前缀保持。新增2AH/2PH、4主+1RK45参考均ground；旧X四raw/五参考仍标旧top460来源。O完整61行×8组合均保存，0–1850km共38有效截面，后续行真实ground后N/A；原600km行保留且通过支持。旧0–550km八组合xyz/p/group差全为0。77项manifest独立核验无漂移。

AH全部实际覆盖及自然端原物理门PASS，L2/PCHIP全程AH−PH最大位置/群程差30.057/7.815mm；AH自然端对本端独立参考约4.261/6.186mm。PH L1/linear s=100km位置1.013233m FAIL，PH原生半径−9.932/−4.422mm及末点path/state严格FAIL保留；其自然端数值门PASS不能覆盖严格接口失败。PH_GRID不解除，不改正确AH事件去拟合PH。本包无新增生产修复。

下一步按[Planner方向](../AH3D-PH474-S5-TOPSUP-20261029-next-direction.md)仅滚动细化WGS84坐标/高度/入口语义；相同ECEF物理位置、波法线和连续Ne/B先对齐，再进入IGRF/空间变B。当前没有下一包实施。完整S5、旧G1/154局部/UNKNOWN/P3.1 0/26及一般截止退化均未因本包解除；旧top460 O支持BLOCKED保留历史。v78原字节见[导航历史](../AH3D-PH474-S5-TOPSUP-20261029/navigation-before-v79/S5-before-v79.md)。

## v78 TOPSUP 实施中快照（历史）

所属[总体v78](../AH3D-PH474-ROADMAP-20261015-plan.md)。2026-10-29状态 `TOPSUP_PLAN_ACCEPT / ENGINEERING_AND_REFERENCE_IN_PROGRESS`。当前唯一细化包[TOPSUP v1](../AH3D-PH474-S5-TOPSUP-20261029-plan.md)，依据[继续实施授权](../AH3D-PH474-S5-TOPSUP-20261029-execution-authorization.md)，工程/参考进行中，尚无新AH/PH轨迹。

只扩同解析T1的顶部支持至600/inner590km，原Ne公式460尺度/真实Ne/B/发射/常数/门槛不变。基本新增O 4+1参考、2AH/2PH；旧X按旧top460来源只读复用，不伪称新grid执行。先bitwise前缀及全部支持消费者自检/参考→工程独立门→原定采集→比较终审。完整O每50km与自然端分别判门，不预设反射，不因终止抹去旧600行。

完整S5/PH_GRID及原局部、UNKNOWN/P3.1边界保持；当前包验收后才滚动细化下一方向。下方v77为旧有限包真实历史。

## v77 常磁场有限包已验收（历史）

所属[总体v77](../AH3D-PH474-ROADMAP-20261015-plan.md)。2026-10-29状态 `CONSTB_OX_LIMITED_DELIVERY_ACCEPT / O_FULL_NATURAL_BLOCKED / PH_STRICT_NATIVE_FAIL`。当前有限包[报告与解决清单](../AH3D-PH474-S5-CONSTB-OX-20261029-report.md)、[独立终审](../AH3D-PH474-S5-CONSTB-OX-20261029-review.md)。18局部、8+2参考、4AH/4PH已完成，生产源码本包未改，124项manifest哈希漂移0。

O只接受0..550km，两端真实事件分别AH physics_domain_rejected与PH -2官方穿透，原600及自然端BLOCKED；X为0..600及AH ground物理PASS。L2/PCHIP AH−PH位置/群程O8.155/4.999mm、X1.319/1.603mm；PH L1/linear两行FAIL、PH X native半径25.856/5.254mm及全path/state FAIL均保留。PH ng字段须乘cos波法线/群方向夹角才对应c/|vg|，该定义已由实际raw确认；不直接将字段沿ds积分。

下一轮只细化[同解析T1扩顶部支持、确认O自然传播](../AH3D-PH474-S5-CONSTB-OX-20261029-next-direction.md)，再按实际验收进入WGS84/IGRF；本轮没有下一包实施。扩域不保证O反射，不重复已接受X全套，只保留必要不变性/同条件控制。完整S5及原PH_GRID/154局部/旧UNKNOWN/P3.1等阻塞不变。

## v76 工程快照（历史）

所属[总体v76](../AH3D-PH474-ROADMAP-20261015-plan.md)。2026-10-29当前状态 `CONSTB_OX_PLAN_ACCEPT / ENGINEERING_IN_PROGRESS`。用户已授权继续实施；[当前唯一细化包v1](../AH3D-PH474-S5-CONSTB-OX-20261029-plan.md)、[独立计划审查](../AH3D-PH474-S5-CONSTB-OX-20261029-plan-review.md)、[授权](../AH3D-PH474-S5-CONSTB-OX-20261029-execution-authorization.md)。同T1真实Ne、40μT常ECEF B，张量O/X与局部物理门→参考资格→工程独立门→有限4AH/4PH→比较终审，正常条件衔接不重复询问。尚无本包AH/PH新轨迹。

AH/PH保留各自plasma/gyro常数；[PH gyro静态物证已落实](../AH3D-PH474-S5-CONSTB-OX-20261029-gyro-static-report.md)。磁化群程用c/|vg|，不套用B=0的1/n；显式检查Domega符号与分支，不把有限安全域结果外推截止退化。旧PH末点、PH_GRID、G1与P3.1阻塞保持。完整S5尚未完成。

## 前置 DIAG3D-B0 历史（v75）

所属总体：[总体计划 v75](../AH3D-PH474-ROADMAP-20261015-plan.md)。2026-10-29 当前状态：`DIAG3D_B0_T1_LIMITED_ACCEPT / G1_BLOCKED`；原两场目标及完整 S5 未完成。

## 当前唯一细化工作包

[AH3D-PH474-S5-DIAG3D-B0-20261029 v1.2](../AH3D-PH474-S5-DIAG3D-B0-20261029-plan.md)，依据[本轮执行授权](../AH3D-PH474-S5-DIAG3D-B0-20261029-execution-authorization.md)，已完成[工程实际ACCEPT](../AH3D-PH474-S5-DIAG3D-B0-20261029-engineering-review.md)及[限定T1最终独立ACCEPT](../AH3D-PH474-S5-DIAG3D-B0-20261029-review.md)。前置 [S4 PHYSALIGN 最终 ACCEPT](../AH3D-PH474-S4-PHYSALIGN-20261029-review.md)有效。详见[总报告与解决清单](../AH3D-PH474-S5-DIAG3D-B0-20261029-report.md)。

原设计含径向G1及多仿射三维T1，B=0、无碰撞、球面，保留真实Ne和物理SI默认。参考8+2已完成；G1/PH loose壳残差1.950862e−10超过原1e−10门，G1未采且继续阻塞。T1已取得两档AH/PH完整raw，AH两重建截面及自然末点物理门PASS，PH注册截面条件门通过，固定共同段偏折约2.516km。两档AH对SI自然端位置约1.72mm/群程3.07mm；L2/PCHIP AH−PH截面位置3.052mm/群程3.549mm。不能把T1接受写成原两场全部完成。

比较分别给出AH−PH、AH−SI、PH−本端条件参考及两参考差，固定13个共同截面和linear/PCHIP均保存。PH内部RK查询不可见，支持仍为 `CONDITIONAL_INTERIOR_SUPPORT`；PH原生地面半径−8.746/−4.090mm及含末点全path/state门FAIL保留，不解除PH_GRID或完整自然末点资格。MATLAB按[上下文恢复补充](../AH3D-PH474-S5-DIAG3D-B0-20261029-matlab-context-recovery-plan.md)完成，累计启动3含原失败1，成功raw2均有真实退出来源。生产/安装包/偏好未改。

## 后续滚动方向

下一主线建议为远离截止/退化的常磁场O/X：先分支身份、局部色散/导数与物理参考，再有限轨迹，待用户确认后细化下一层，本轮不实施。见[Planner反馈](../AH3D-PH474-S5-DIAG3D-B0-20261029-next-direction.md)。G1重复控制、PH末点及独立机制份额在影响后续验收或明确修复决策时恢复处理，不以它们未穷尽为由阻塞无依赖的主线。原154局部阻塞不外推解除。WGS84/IGRF继续后续大纲。

有限前缀、人工截停和上边界退出不能冒充完整自然反射。P3.1 原 26 例与 0/26、154 局部阻塞、PH 原末点壳面 FAIL、旧 UNKNOWN 和严格完整资格缺口保持。最终结果供 [S6](./S6-plan.md)汇总。旧页面原字节见[导航历史](../AH3D-PH474-S5-DIAG3D-B0-20261029/navigation-before-v74/S5-before-v74.md)。

