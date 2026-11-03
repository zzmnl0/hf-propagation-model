# AH3D 生产与科研统一 L2：迁移结果

Task：`AH3D-PH474-L2-DEFAULT-20261103`；2026-11-03。用户直接授权“将生产和科研等均调整为使用L2，后续也均采用L2”。本报告记录实施与验证；实际独立决定见[终审](./AH3D-PH474-L2-DEFAULT-20261103-review.md)。

## 当前行为

- `config.py`以`AH3D_L2_PROFILE/SETTINGS`保存统一L2；默认ECEF模型、fan、P2P积分及最终回算、main2 CLI均采用L2。正式入口拒绝显式低档及非canonical L2设置，底层显式诊断档保留。
- L2为rtol=1e-10、位置atol=1e-8km、归一化波矢atol_p=1e-12、步长.002–.5km，普通空间FD基准.03125km。细结构沿既有保护继续细化FD，实际值记录，不伪称逐参数等于原冻结场景。
- 生产积分状态存物理k，按每条轨迹频率转换`atol_k=k0*atol_p`并记录实际参数；旧显式atol_k保留原单位，同时指定两单位明确报错。
- P2P搜索所得轨迹已是L2时沿原流程复用；真正发生的最终回算也只能使用L2，不为同档重复积分。P2P收敛门、几何、介质、事件和物理阈值未放宽。
- 未来生产和科研结果统一L2。历史低档只用于明确诊断和版本复核；原SIACCESS绑定旧档的测试原样逻辑归档，新L2回归承接当前路由。正式轻量runner现为七项，未以glob运行历史测试。

实现及后续入口见[脚本索引](../AH3D-SCRIPTS.md)、[当前测试索引](../tests/README.md)和[Worker报告](./AH3D-PH474-L2-DEFAULT-20261103/worker-report.md)。

## 实际验证

- 仅选择正式runner中的`_verify_l2_defaults.py`，实际PASS。8条各.25km短轨迹：8/12MHz真空O/X对解析直线及群程回归，另平滑磁化O/X与同进程显式冻结L2设置逐值一致。实际DOP853参数确认atol_k/k0=1e-12；路由替身与真实积分分开记录。
- 当前`_verify_ah3d_physalign.py`单独PASS，额外2条各.05km短轨迹，检查SI/显式系数及兼容控制。实际总计10条短轨迹，不把计划建议8条伪报成实际总数。
- CLI帮助/选择、fan/P2P/最终回算路由、正式低档拒绝、非法设置、FD保护和兼容单位已针对性检查。两受控进程均正常退出，单线程，未触发资源停止；日志与结果在本Task新目录。
- 完整短轨PASS后，唯一核心补改为缺省feature_gates显示标签research→L2，对应`settings_checks()`定点PASS。完整8轨运行早于该标签补改；实际积分/参数实现没有再改。其后的runner计数文档更正不重采轨迹。
- 静态确认核心没有accepted-step数量上限；正式fan已有8000km路径上限，不会因.5km步长触发旧4000步截断。P2P max_nfev是优化求值次数，未改；长路径性能仍未实测。

这证明配置/入口迁移和有界回归，不等于新一般三维、完整P2P、长自然端或散射功率验收。未启动PH/MATLAB/IRI/IGRF/REA，未重跑旧26例或冻结PH数据。原FAIL/UNKNOWN及各正式矩阵状态保持。

## 原版本保留与身份

HEAD仍为`4be052c0c69f31a631beb8ed823c13d539036274`，以工作区文件哈希区分本次版本，不能仅用HEAD表示新代码。原15项tracked修改未由HEAD覆盖。

[before映射](./AH3D-PH474-L2-DEFAULT-20261103/before-worker/manifest.json)保留10个既有文件的原字节，[after清单](./AH3D-PH474-L2-DEFAULT-20261103/after-worker-manifest.json)绑定当前文件，[差异补丁](./AH3D-PH474-L2-DEFAULT-20261103/worker-diff.patch)记录改动。[历史版本核验](./AH3D-PH474-L2-DEFAULT-20261103/historical-version-resolution.json)确认：

- 前一WGS科学528叶中522项当前仍匹配，6项已授权源文件升级可由before原字节精确复核，无未解释漂移。
- AH清理13叶中11项当前匹配，2项入口索引升级可由before快照精确复核，无未解释漂移。
- 旧科学数据、冻结脚本和manifest保持原件，不宣称旧清单对迁移后的活动源路径仍零漂移，不回填旧通过数。

总体/S5/S6的v100原导航已保存，后续导航仅登记L2政策与本次真实审查。当前新绑定与行政收尾见[候选清单](./AH3D-PH474-L2-DEFAULT-20261103/candidate-manifest.json)和[收尾](./AH3D-PH474-L2-DEFAULT-20261103/closeout.json)。

## 后续

后续生产及科研AH轨迹统一使用L2；不再以低档生成正式结果。计算成本预计增加，尤其fan/P2P；本轮未作完整性能承诺。下一科学方向仍为S6集中汇总已有误差定位、已修/待修/未确认项及反射链路适用域，不自动开展新科学包。
