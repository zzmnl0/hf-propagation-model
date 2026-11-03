# AH3D 生产与科研 L2 统一迁移计划 v1

Task：`AH3D-PH474-L2-DEFAULT-20261103`；日期：2026-11-03；项目：`D:/code11/IRI01/IRI03/TID/code`。由实际Planner `/root/l2_default_plan`制定，Root保存。用户已直接授权迁移，见[真实授权](./AH3D-PH474-L2-DEFAULT-20261103-authorization.md)，不再申请相同范围批准。

## 1. 目标与边界

将已验收L2配置迁入统一配置，覆盖默认ECEF模型、fan、P2P积分、最终回算及活跃科研CLI。正式入口不再产生低档输出；旧档只供底层显式诊断、负控制和历史复核。沿用SI普通模式默认和显式系数接口。

L2不改变地球几何、真实Ne/B、支持域、地面/顶部、返回分类、P2P收敛门或物理阈值。本次精度迁移不自动取得一般场景、P2P完整性或散射功率资格。当前HEAD `4be052c0c69f31a631beb8ed823c13d539036274`，迁移前15项tracked既有修改必须保留。

## 2. 已核实断点与范围

| 文件 | 当前行为 | 迁移 |
|---|---|---|
| `config.py` | 无稳定L2，总档research | 唯一canonical L2及普通默认 |
| `models/ray_tracer_ah_3d_ecef.py` | 默认reference/FD .5，积分状态k | L2、FD .03125、p容差按k0换算 |
| `models/unified_hybrid_pipeline.py` | fan smoke、factory FD .5、P2P reference | 统一L2，公共正式入口拒绝低档 |
| `models/point_to_point_ah_3d_ecef.py` | 默认reference | 省略时随模型L2，显式底层诊断保留 |
| `main2.py` | CLI只接受旧档、摘要fallback reference | 选择/默认/透传/摘要统一L2 |
| `main.py` | 已走统一入口 | 用路由测试确认，无预定编辑 |
| `analysis/` | 无直接AH追迹/P2P或旧档调用 | 无预定编辑 |

`solve_profile`和`final_retrace_profile`在当前ECEF链均为轨迹积分档。P2P的tol_km、Jacobian、多尺度搜索、去重及候选完整性是其他语义，不能因L2自动宣称重新验收。

## 3. canonical L2与单位

来源：已接受REFLECTWGS的`wgs_ah.py`和`AH/O_L2/request.json`。

| 参数 | 值 |
|---|---:|
| rtol | 1e-10 |
| atol_r_km | 1e-8 |
| atol_p | 1e-12 |
| ds_min_km | .002 |
| ds_max_km | .5 |
| 普通模式名义spatial_fd_base_km | .03125 |

生产状态保存物理k，实际应为`atol_k=model.k0*atol_p`。不得直接用atol_k=1e-12或固定12MHz换算值。配置存atol_p，积分入口按频率换算；原显式atol_k语义保留，同时提供两单位应明确拒绝。新容差须有限正值，并保存实际消费参数或等效运行证据。

factory保留`min(base_fd, minimum_structure_scale_km/60)`，仅把普通base从.5改为.03125；细结构可更小，记录实际值。FD scales/stability、投影1e-10、事件定位1e-9km、群时状态容差及投影开启行为沿用。SI系数80.61638587963631已正确，不重复修；显式及旧兼容模式保留追溯语义，不当作新正式L2物理验收。

## 4. 实施顺序

1. 核对Git身份与实际修改列表，保存所有将改文件的当前工作区原字节、原路径、长度和SHA；包括原未dirty的main2/P2P。不能以HEAD内容覆盖既有用户修改。
2. 迁入配置、直接构造默认、普通FD及容差转换；省略profile最终落到统一L2。
3. 贯通fan、P2P积分与最终回算、main2；公共生产/科研明确拒绝显式旧档，不能静默降档。底层保留显式旧档及自定义冻结profile作为诊断。
4. 更新README、AH3D-SCRIPTS、tests/README及总体/S5/S6当前导航，持久记录未来统一L2、成本及验收边界。当前回归清单按实际改动同步。
5. 实际有界验证后完成独立Reviewer审查；只以真实结果接受。

## 5. 验证及验收

新增当前回归`tests/_verify_l2_defaults.py`，用指定CPU Python、单线程、受控超时；输出进入本Task新目录。测试目标：

- 配置、直接构造、factory、fan、P2P及最终回算实际L2；CLI帮助、默认摘要一致。
- 公共低档请求被拒绝，未知profile及非法新容差不静默回退；底层显式旧atol_k语义保留。
- 至少两个频率核对atol_k/k0=1e-12及rtol、位置容差、步长、名义和细化FD。
- 实际短程真空O/X对解析直线和群程回归；平滑磁化介质O/X的canonical L2与同进程显式重建冻结L2设置一致。后者是迁移等价，不冒充新独立物理参考。
- 公共分发/种子/P2P参数透传可使用可观察替身，实际积分必须另行真实执行，二者分开记。建议总计最多8条、各≤10km短轨迹，不运行全场P2P、完整主场景或旧26例。
- 合理有界进程，失败定位后只重跑相关项；不采冻结PH或长自然轨迹。

现有`_verify_ah3d_physalign.py`为当前通用回归，其中.5km默认断言应随迁移更新并保留before快照；继续检查兼容及显式系数优先级。

`_verify_siaccess_request.py`绑定旧任务帮助模块、旧清单/构造额度及smoke/research/reference路由，标为历史冻结入口，由新L2测试承接当前透传；不修改旧脚本制造旧验收继续通过。

验收：新增当前回归实际通过；没有未解释的正式低档默认/绕行；不放宽原门，不丢用户修改；区分配置与短程验证通过、未重跑完整科学范围；实际Reviewer独立作决定。

## 6. 历史版本与交付

历史科学清单绑定config、README及生产源码，清理清单绑定AH3D-SCRIPTS。授权迁移导致这些活动路径的预期版本差异：

- 旧科学数据、manifest、报告、归档脚本原件保留。
- 提供原路径/哈希到before快照的映射，按旧代码版本复核。
- 不宣称旧manifest在迁移后活动路径仍零漂移，不将预期升级记作旧科学失败，不回填旧通过数。
- 新结果绑定迁移后代码及测试证据。

交付计划、授权、before快照/映射、Worker报告、实际测试记录、独立review及收尾。总体记录L2政策和新代码身份。无需再询问正常适配与有界修复；若需改变物理场景、放宽门、改旧证据或启动未授权原生程序，停止对应分支并报告。
