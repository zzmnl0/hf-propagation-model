# S2 判别关键物理条件：阶段入口

- 所属总体：[总体主计划v26](../AH3D-PH474-ROADMAP-20261015-plan.md)；本入口版本v21；状态：`S2B_ENTRANCE_PREFIX_DELIVERY_COMPLETE / PH_ENDPOINT_EVIDENCE_PACKAGE_ACCEPTED / NATURAL_ENDPOINT_GATE_OPEN / NATEND_R1_E_REPLAN / NATEND_R2_E_REWORK / NATEND_R3_E_REPLAN / NATEND_R4_E0_REWORK / NEXT_EXECUTION_NOT_AUTHORIZED / RETRY_NOT_AUTHORIZED / NOT_GENERAL_EQUIVALENCE_ACCEPTED`。
- 最近实施：[NATEND-R4来源准入E0最终v1](../AH3D-PH474-S3-NATEND-R4-20261017-plan.md)已获[明确执行授权](../AH3D-PH474-S3-NATEND-R4-20261017-execution-authorization.md)，但唯一尝试为实际独立 **REWORK**；见[执行报告](../AH3D-PH474-S3-NATEND-R4-20261017-execution-report.md)、[实际审查](../AH3D-PH474-S3-NATEND-R4-20261017-execution-review.md)及[外层回执](../AH3D-PH474-S3-NATEND-R4-20261017-execution-receipt.json)。ParseFile通过，初始悬挂子入Job为1/1、Resume后计数门首失；触发计数未落盘，清理复用同门，最终Job活动数及自动成员身份UNKNOWN。已知子exit99、外层工具exit1/no session不等于全Job退出或E0通过。子端双向核验/实际错对象及模拟控制未取得完成证据，显式Python/科学调用0不能认证未知自动成员的类型或数量。

本次已消耗唯一尝试，余窗不重开；完整E、D2及R/S/C仍未放行。当前不新细化返工包。下一方向仅为计数先保存后判断、清理独立记录最终活动数与有界成员身份取证，另行规划和许可；不预设辅助成员豁免或放宽上限，不据工程失败新增AH/PH物理错误。

更新前总体v25及S2/S3/S4原字节见[执行历史](../AH3D-PH474-S3-NATEND-R4-20261017/execution-history/)，[执行导航版本映射](../AH3D-PH474-S3-NATEND-R4-20261017-execution-navigation-version-map.json)仅按原路径及期望SHA将本任务v25引用解至同哈希历史；旧v24及更早映射、冻结数学及原科学FAIL/UNKNOWN保持。历史相对链接依原目录语境解释。
- 更新前总体v24及S2/S3/S4原字节见[R4历史](../AH3D-PH474-S3-NATEND-R4-20261017/history/)，[R4显式版本映射](../AH3D-PH474-S3-NATEND-R4-20261017-navigation-version-map.json)仅将R3执行回执27中的4导航及旧执行映射/收尾审查的4项current-v24引用解至同哈希历史；其余23项原路径、v23及更早解析不变，不改旧证据。历史相对链接依原目录语境解释。
- 前序实施（R3历史时点）为NATEND-R3最终v1，[实际独立REPLAN](../AH3D-PH474-S3-NATEND-R3-20261017-execution-review.md)，见[报告](../AH3D-PH474-S3-NATEND-R3-20261017-execution-report.md)。PID0控制NULL/87通过，来源父唯一0x1000查询NULL/5；具体策略及R2旧原因仍未确认。Job/controller/Python/D2及新增科学均0，R/S/C未开，R3停止时尚未新细化工作包，当时仅建议重规划可核验的启动来源链，不提权、不跳过身份门或使用余时重试。实施前v23四导航见[执行历史](../AH3D-PH474-S3-NATEND-R3-20261017/execution-history/)，[执行版本映射](../AH3D-PH474-S3-NATEND-R3-20261017-execution-navigation-version-map.json)只解planning30的4/26、prelaunch68的4/64及旧映射4项current-v23；其余旧来源/数值不改。
- 前序实施（R2历史时点）：NATEND-R2最终v1已获明确批准，唯一E包装于2026-10-17T11:29:00Z因父进程OpenProcess返回0硬停止，见[执行报告](../AH3D-PH474-S3-NATEND-R2-20261017-execution-report.md)、[实际独立REWORK](../AH3D-PH474-S3-NATEND-R2-20261017-execution-review.md)及[外层回执](../AH3D-PH474-S3-NATEND-R2-20261017-execution-receipt.json)。Job/controller、Python/guard/D2及新增数学均0；R/S/C和真实科学活动仍0，不重新开放S2A/S2B采集。下一方向仅最小只读访问修正及LastError留证；实际因果及最小query成功仍UNKNOWN，重试须另行补充计划和许可，不用余时自动重开。
- 更新前本入口v16及总体v21四导航见[执行历史](../AH3D-PH474-S3-NATEND-R2-20261017/execution-history/)；[执行导航映射](../AH3D-PH474-S3-NATEND-R2-20261017-execution-navigation-version-map.json)仅重定位R2 planning25中的4导航及旧映射4个current绑定，其余21项和全部E/数值来源保持原路径。更早[v15副本](../AH3D-PH474-S3-NATEND-R2-20261017/history/S2-plan-before-natend-r2.md)及[R2原映射](../AH3D-PH474-S3-NATEND-R2-20261017-navigation-version-map.json)不改；下文旧记录保留历史时点语义。
- PHEND R1的[实际独立ACCEPT](../AH3D-PH474-S3-PHEND-R1-20261017-execution-review.md)仅接受终点语义/重放证据；两个原接口FAIL继续保留。本轮将未关闭的终点科学门连接到[S3 NATEND独立完整自然终点基准计划](../AH3D-PH474-S3-NATEND-20261017-plan.md)，不重新开放S2A/S2B采集。更新前[v13原字节副本](../AH3D-PH474-S3-NATEND-20261017/history/S2-plan-before-natend.md) SHA `EBB8459A2C727A0E1837884098B32280FB4E25DE5AA489B28CD59C07385A7967`；[新版本映射](../AH3D-PH474-S3-NATEND-20261017-navigation-version-map.json)保存旧PHEND规划绑定。更早v12及FINISH绑定见[旧映射](../AH3D-PH474-S3-PHEND-20261017-navigation-version-map.json)，旧文件不改。
- 最近详细任务：[FINISH恢复v1](../AH3D-PH474-S2B-FINISH-20261017-plan.md)依[直接执行授权](../AH3D-PH474-S2B-FINISH-20261017-authorization.md)及[独立E ACCEPT](../AH3D-PH474-S2B-FINISH-20261017-engineering-review.md)完成唯一剩余PH。10项诊断齐备、8组分解完整，S/C实际exit0，历史AH4/PH8；[执行报告](../AH3D-PH474-S2B-FINISH-20261017-execution-report.md)与[最终实际审查](../AH3D-PH474-S2B-FINISH-20261017-execution-review.md)独立绑定，不追认原全返回门或一般S2B通过。下文为保留的历史步骤。
- [v11入口原字节历史副本](../AH3D-PH474-S2B-FINISH-20261017/history/S2-plan-v11.md)，SHA `0A8024E3340CA6597021F655A035FB07440DCB397B68C38E63805B3F78A2F7AF`；旧审查不自动接受本v12。
- [v10入口原字节历史副本](../AH3D-PH474-S2B-PREFIX-20261016/history/S2-plan-v10.md)，SHA `53BD7B439633BAD7A2556DAA14D3BB7679CA7470526DD05599CFAAB377B70010`；旧PLAN审查仍绑定历史v10，不自动接受本v11。
- [v9入口原字节历史副本](../AH3D-PH474-S2B-PREFIX-20261016/history/S2-plan-v9.md)，SHA `807831C7C069EF6114CF0B61129899B16FBEEF34D243AC0FCB5A716BA2EEE2FB`；前序实际审查仍绑定历史v9，不自动接受本v10。历史相对链接依原目录解释。
- [v8入口原字节历史副本](../AH3D-PH474-S2B-REPAIR-20261016/history/S2-plan-v8.md)，实际SHA `9F9CBEB77763F52D3598A815AB804E14474CEDE0E3AF9F626DB8933BEDA4BC9A`。旧恢复execution-review文件实际为SHA表，所列S2哈希不同于本次快照实字节；不覆盖旧文件，不以旧表追认本v9，新审查独立绑定当前版本。
- [v7入口原字节历史副本](../AH3D-PH474-S2B-RESUME-20261016/history/S2-plan-v7.md)，SHA `AF2BF5492B95E07738D5EB94E3062FB8BB52BA9EDCC74AA33F1A52BFD7067039`；旧审查仍绑定历史v7，不自动接受本v8。
- [v6入口原字节历史副本](../AH3D-PH474-S2B-20261015/history/S2-plan-v6.md)，SHA `D6522782778D0AA36312E8F8F29A43BED6BC9269FEDAE35AD0ED1F94129F5106`；历史相对链接依原目录解释，旧计划审查不自动接受本v7。
- R4最终独立结论：`ACCEPT`，仅接受诊断交付及冻结常介质有限状态的局部一致性；不是无条件Gate A、完整S2、AH–PH等价或连续域验证。
- [v5入口原字节历史副本](../AH3D-PH474-S2B-20261015/history/S2-plan-v5.md)，SHA `4A213819EEE5137D016FF5FDE6C9A269AB399031670FAB5E551B715DABCE9322`；历史相对链接依原目录解释。R4文档审查仍绑定v5；本v6与S2B计划由[新独立计划审查](../AH3D-PH474-S2B-20261015-review.md)另行评定。
- [v4入口原字节历史副本](../AH3D-PH474-S2A-ENTRY-R4-20261015/history/S2-plan-v4.md)，SHA `A308E3765CA4616E148559048E394777E194C272C71A8BC9FCD80EDCD51A8798`；历史相对链接依原目录解释，旧审查不自动接受本v5。
- [v3入口原字节历史副本](../AH3D-PH474-S2A-ENTRY-R3-20261015/history/S2-plan-v3.md)，SHA `407D89063EB98B8E5AB12BEAC7F4D4C0E172D01700156A1D34DDF8275EF701A3`；历史相对链接依原目录解释，旧审查不自动接受本v4。
- [v2入口原字节历史副本](../AH3D-PH474-S2A-ENTRY-R2-20261015/history/S2-plan-v2.md)，SHA `72A357A2FAA2262D14F6EAB5CDFA685FFEC60D03B5F21C5438404CCB9E658F53`；历史相对链接依原目录解释，旧审查不自动接受v3。
- [v1入口原字节历史副本](../AH3D-PH474-S2A-ENTRY-R1-20261015/history/S2-plan-v1.md)，SHA `1EB11980F0F0DF639C5773EF8BC5938047E9E8D00A7626AAC89C925224FD3521`，来源为本文件更新前版本；其中相对链接依原所在目录解释。旧审查绑定v1，不自动认证v2。

## 当前完成层级

[R4最终实际审查](../AH3D-PH474-S2A-ENTRY-R4-20261015-execution-review.md)已接受限定结果：32新mode/64方向/128新状态，加R3六个合格完整mode/24历史状态，共152唯一状态与76旋转检查，保存判断全部PASS。旧I003_O两条部分状态保留、排除于主分母，本次已完整重新采集该mode。154原阻塞、原P3.1契约与0/26均不改变，不写成192项或一般三维通过。

以下为已保留的历史步骤，不代表需要重做：

现有参考方向、条件预测及R1控制可按哈希/接口兼容复用。随后旧S2A实施已补emit并曾获E ACCEPT，但[最终实际审查](../AH3D-PH474-S2A-SCI-20261015-execution-review.md)为REWORK：首模式前超时，0科学观测；另确认真实状态封装与collector不兼容。该schema缺陷未实际触发，具体导入耗时原因仍未知。历史E接受和失败证据均保留，不能作为新入口已就绪的依据。

[ENTRY-R1实际审查](../AH3D-PH474-S2A-ENTRY-R1-20261015-execution-review.md)随后确认schema/完整collector/C工程进展，但两次就绪探针均在来源alias审计拒绝，最终REWORK，无S窗口。下一步保留完整工程证据，集中补齐固定来源映射和全模块origin/hash/match观测；Probe1缺失origin仍UNKNOWN。

[ENTRY-R2实际审查](../AH3D-PH474-S2A-ENTRY-R2-20261015-execution-review.md)为REPLAN：完整844记录中16项拒绝，未启动S。[R3实际E审查](../AH3D-PH474-S2A-ENTRY-R3-20261015-engineering-review.md)随后ACCEPT：有限九成员/七对象规则、154差量控制、唯一probe 11.882229秒/844全匹配及42成员回执成立；旧UNKNOWN不回填。

[R3报告](../AH3D-PH474-S2A-ENTRY-R3-20261015-report.md)与[最终实际审查](../AH3D-PH474-S2A-ENTRY-R3-20261015-execution-review.md)记录一次真实S在I003_O部分状态因新增输出资源门停止。完整mode为R007/I001/I002各O/X；另I003_O rotation0两壳比较已保存。4正12负逐项覆盖和12旋转记录保留，即时checker嵌套汇总false的语义不改写。正式科学汇总缺6文件，89成员行政封存不是科学接受，也没有证实AH物理缺陷。

[R4报告](../AH3D-PH474-S2A-ENTRY-R4-20261015-report.md)承接上述证据，完成七流无损存储、真实资源基线/退出观测及独立审查后的完整模式复用。首次E闭包审查REWORK及旧42成员保留；修正参考流/调用审计关联后，[43成员实际E复审](../AH3D-PH474-S2A-ENTRY-R4-20261015-engineering-review.md)ACCEPT，再完成唯一S/C。37科学回执及1779清单成员、七个实际压缩流独立解压哈希均核验通过。

## S2A：已完成结果

R4已解决原日志资源阻塞并完成剩余限定采集，不再把已接受的证据工程列为下一轮主线。实际S/C子进程及外层工具均exit0，846来源全匹配、线程1，参考调用128/256/384/12288，未超原门槛。E生产导入/核心/参考数学为0，C只处理已有数据；未改生产AH、PH或数学参考。

局部比较以80位Booker隐式微分中心参考为目标，生产为复数中心差分；指标为归一化后的绝对差/范数，不是相对误差。三h、80/120位FD另作稳定性与分支检查。常数差值实际为reference80中心与cached条件预测的差，非PH实测、非孤立常数贡献。冻结`scientific_acceptance=false/PENDING`字段保留，后出的独立审查给出有限范围接受，不改写历史。

有限连续介质/初态可比性与残差识别方向随后已单独批准并部分实施，见下节。当前结果仍不支持仅为拟合PH条件预测而直接修AH；也不能自动把有限点或常介质单例外推成一般连续轨迹域，不展开S3–S6。

## S2B：历史恢复与字段适配阻塞

目标是受支持条件下 AH、PH 原生轨迹及独立参考的受控比较，确认局部结果如何传到轨迹、事件和自然终止。连续介质及初态的可比性必须有证据；局部有限点不证明连续场。

用户已明确批准[S2B入口v1](../AH3D-PH474-S2B-20261015-plan.md)的两种B=0介质、AH4/PH6上限及E60/S40/C15预算。实际E-B1同窗返工后[独立复审ACCEPT](../AH3D-PH474-S2B-20261015-engineering-rereview.md)，182行参考及33控制复用，不改生产/参考算法。

实际S仅尝试AH1/PH1：C0/L1 AH完整60km返回，13截面linear/PCHIP对本端独立参考均保存PASS；MATLAB随后在启动阶段报`failed to load settings errors_warnings plugin`并exit1，无PH原始MAT。R1/R2及其余6项全部NOT_RUN，不重试、不换参数、不补参考或残差分解。S controller实际exit2，唯一C只封存已有证据亦exit2，23保存成员中缺`native_decomposition.json`明确UNKNOWN/NOT_RUN，527成员最终清单保留原始外层日志/退出。

以上是2026-10-15历史PARTIAL及启动阻塞，原证据不改写。2026-10-16 [沙箱外启动实际PASS](../AH3D-PH474-MATLAB-RETRY2-20261016-report.md)后，用户直接要求按指定PH目录执行未完成部分；[新恢复E独立ACCEPT](../AH3D-PH474-S2B-RESUME-20261016-engineering-review.md)，182行参考和C0/L1 AH原字节复用。

新S唯一一次PH实际exit0、完整返回248节点；不变analysis.py读取ray_path_data.cumulative_absorption发生KeyError。实际path含absorption而无cumulative_absorption，仅确认字段契约不匹配，不猜同义或单位，也不认定PH/AH物理错误。没有timing/initial/comparison落盘，后续8项依赖均NOT_RUN。新尝试AH0/PH1、历史AH1/PH2，完整逻辑结果仍仅复用AH1/10。新S/C外层exit2，纯封存partial，不重试、不补算参考/比较或分解。

PH公开raytrace_3d_sp.m:241也将该path字段写作cumulative_absorption，实际MEX返回名称不同；因此不能把问题归为Worker凭空使用字段，也不能仅据absorption名称确认别名语义。新71成员manifest及27保存成员receipt、旧527成员均已复核，最终结论由本轮独立审查绑定。

以上是v8时点状态；该诊断适配返工现已在独立新任务中实施，结果见下节，不追认旧运行成功。

固定截面残差不是连续前缀误差上界，PH公开输出不提供内部原生Ne/梯度或D_omega。PH无已确认固定几何程硬停接口，完整返回与用于比较的内点前缀分别记录；不把资源截断或输出标签单独当自然事件等价。详细门、参考和停止条件以当前计划为准，旧S2 v2的B部分仍仅是历史候选。

当前D_omega符号条件风险属于AH弧长参数化，未确认实际触发；PH对应一般代数未完全确认，不命名为已确认PH守卫bug。B=0解析负号不替代一般磁化数值筛查。仅为未来S3保留可识别性需求，不在本轮细化或承诺投影/积分/事件唯一份额。

## S2B：本轮适配修复后仍按接口门停止

[本轮E实际独立ACCEPT](../AH3D-PH474-S2B-REPAIR-20261016-engineering-review.md)接受75成员、18控制及旧527+71原字节来源核验。原生path/state均为absorption，原全部累计量、初态/时刻及参考阈值保持。E无生产/MATLAB/真实PH校准/参考数值调用。

唯一S先复用已保存C0/L1 PH raw完成原时刻/初态/13截面两重建门，再完成R1/R2及C0/L2两端、G1/L1两端的六次新调用。六新child均exit0，C0及restart保存PASS，G1/L1 AH保存PASS；新AH2/PH4、历史AH3/PH6。G1/L1 PH虽完整返回1174节点，但全数组位置对应最大差3.4357485875e-05km>1e-8km，原门保持STATE_TIMING_UNRESOLVED。2项G1/L2依赖NOT_RUN，7/10比较完成、8/10原始返回，非S2B完成。

[独立接口诊断](../AH3D-PH474-S2B-REPAIR-20261016-g1-interface-review.md)确认只在最后ground节点超门，几何程1731.814km、path.height=0而state半径高度约-3.44cm；角度一致。其余节点及0-60km前缀位置对应观察在门内，不替代缺失的G1/PH科学比较。底层事件处理为推测，原门无末点豁免，独立决定REPLAN。不能删末点或投影state把原v1改判PASS。

S实际08:20:18.1336538–08:21:12.1680345Z退出2；C实际08:27:23.3720635–08:27:24.2917870Z退出2，C不补目标/参考/比较。173 manifest与100/101科学成员哈希核验通过，唯一缺native_decomposition。原182行参考未重算，仅原许可restart新增10行常量参考、仿射积分0；未用窗口/额度不恢复。

用户已选择只规划前缀/终点分门的补充方向，见[真实选择记录](../AH3D-PH474-S2B-PREFIX-20261016-direction.md)与[契约研判](../AH3D-PH474-S2B-REPAIR-20261016-contract-triage.md)。前序独立实际REPLAN之后，本轮仅细化该工作包并独立PLAN审查；方向选择不授权新契约或续跑。

## S2B：分门补充计划与实际partial

[PREFIX补充计划](../AH3D-PH474-S2B-PREFIX-20261016-plan.md)是唯一详细下一步。拟复用8raw、7比较、192参考行，零新参考，仅新增G1/L2 AH1/PH1，历史累计最多AH4/PH7；拟E30/S20/C10新窗口，独立E实际ACCEPT后才条件进入S。计划接受与用户执行批准分开，不恢复旧窗口或剩余次数。

固定主查询0-60km、restart旧相对查询；linear/PCHIP用完整原数组重建，接口验证覆盖全部依赖包络。必须证明整个扩展包络及AH有限差分域的介质支持；PH事件须在包络及至少一个内域缓冲节点之后。证据不足停依赖，不通过截短数组、缩短截面或新增参考求积取得通过。

原全返回结构、累计量及位置门保留。仅原生完整单hop、已知末ground径向异常且包络/缓冲外等全部严格条件成立时，允许前缀另行判定；全返回仍FAIL、底层机制UNKNOWN、自然终点等价NOT_ESTABLISHED。10逻辑案例的完整诊断交付不等于残差全部PASS，也不追认原S2B最终v1通过。

以上为v10规划时内容，现已明确批准并实施，不把后出的授权回填进冻结计划。新E52+5控制、178成员独立ACCEPT；唯一S复用7比较后，新完成G1/L1 PH前缀门/13截面比较与G1/L2 AH比较，9项均保存本端有限参考PASS。G1/L1 PH包络到62.8014770835km、J41/buffer42/首事件1173，原末点全返回FAIL仍保留，不认证自然终点等价。

唯一G1/L2 PH在已启动MATLAB的s2b_export第20行raytrace_3d_sp调用处报“内存不足”，child exit1且无raw。采样未超4GiB树RSS或512MiB输出帽，无控制器终止动作；消息不证明系统总内存耗尽，底层原因未知。新AH1/PH1、历史AH4/PH7已用完，不重试。S09:15:55.8191107–09:16:19.6046109Z、C09:18:18.3071746–09:18:20.8653434Z均outer exit2。

C纯封存226成员manifest及163/168科学成员，明确缺G1/L2 PH raw/initial/prefix/comparison与native_decomposition五项；192参考行零重算，停止后科学工作0。输出已冻结，生产/PH包未改。最终审查后仅建议限定定位本次原生失败，待用户确认范围/预算再细化，不进入S3-S6。

## 转出条件

当前S2证据仍以[FINISH实际报告](../AH3D-PH474-S2B-FINISH-20261017-execution-report.md)和[最终审查](../AH3D-PH474-S2B-FINISH-20261017-execution-review.md)为准：第10项原始返回/前缀/比较及8组条件代数分解已保存，九旧比较原字节复用；140清单成员、232来源引用受哈希约束。G1/L2 PH2344节点，NRT=1的唯一调用成功，不能据此认定旧OOM根因。原L1约3.44cm及新L2约8.317mm末地面节点均为全返回FAIL，前缀PASS不认证自然终点等价。S3A及SIRESID有限离线归因随后已接受，当前唯一下一步是[S3 PHEND计划](../AH3D-PH474-S3-PHEND-20261017-plan.md)，优先已保存末段与静态证据；不新增S2调用，执行仍待批准。若后续形成AH独立物理缺陷依据，仍可选[S4](./S4-plan.md)最小修复，不要求先完成PH全域等价。

本入口v13只更新后续导航，不新增科学运行。原失败、UNKNOWN及154阻塞保留；D_omega符号风险、PH一般代数、截止/退化、投影/积分/事件独立贡献、一般三维/O-X连续场、完整自然终点及WGS84/IGRF仍待验证。未修改生产AH或PH包、未重跑正式P3.1（仍0/26）；不回写旧false/PENDING，不将本端有限PASS写成AH-PH直接差消失。

