# AH3D L2 默认迁移实施记录

任务 `AH3D-PH474-L2-DEFAULT-20261103`，计划 v1。用户授权原文和范围见[授权记录](../AH3D-PH474-L2-DEFAULT-20261103-authorization.md)。基线和当前 Git HEAD 均为 `4be052c0c69f31a631beb8ed823c13d539036274`；开始时工作区已有 15 项 tracked 修改。实施前对拟修改文件保存了工作区原字节及哈希，没有以 HEAD 覆盖这些修改。原文件映射见[before 清单](before-worker/manifest.json)，当前文件长度和 SHA-256 见[after 清单](after-worker-manifest.json)，本工作包的文件级差异见[差异补丁](worker-diff.patch)。before 中现存文件的快照长度和 SHA 已逐一核验通过。

## 实施内容

- `config.py`、ECEF 模型、统一入口、ECEF P2P 和 `main2.py`：正式默认及普通接口采用 L2，`rtol=1e-10`、位置 `atol=1e-8 km`、无量纲波矢 `atol_p=1e-12`、步长 `0.002–0.5 km`；普通空间差分基准为 `0.03125 km`，细结构允许继续缩小。
- 积分状态保存物理波矢，故每条轨迹按其 `model.k0` 将 `atol_p` 转为 `atol_k=model.k0*atol_p`；返回元数据记录实际消费参数。底层显式旧 `atol_k` 档用于历史诊断；同一档不能同时指定 `atol_k` 和 `atol_p`。
- 正式 fan/P2P/CLI 拒绝旧档；底层显式历史档及 PH 径向兼容配置保留。P2P 优化及最终档同为 L2 时，复用已收敛轨迹并记录 `same_profile_solution_reused`，不额外重复回算。
- 更新当前回归清单、新增 L2 回归、调整已有物理系数路由测试，并更新 `README.md`、`AH3D-SCRIPTS.md`、`tests/README.md`。统一入口 `feature_gates` 的缺省档标签也与 L2 一致。
- 文档收尾将 `AH3D-SCRIPTS.md` 的正式 runner 项数更正为七项，并区分历史清理当时未迁移默认与本次另获授权的 L2 迁移；此后仅刷新清单和差异记录，未重跑测试或修改核心代码。

## 实际验证

以下命令运行于 `D:\code11\IRI01\IRI03\TID\code`，使用 `C:\Users\12238\.conda\envs\pytorch_cpu\python.exe`。回执中记录了线程环境、退出码和采样进程树资源；使用新输出目录重现时须避免覆盖本次产物。

1. `python -B tasks/AH3D-PH474-L2-DEFAULT-20261103/bounded_checks.py l2`：**PASS**，正式 runner 仅选择 `_verify_l2_defaults.py`。见 [l2 回执](l2-receipt.json)、[runner 摘要](runner-l2/summary.json)和[数值明细](runner-l2/_verify_l2_defaults/l2-defaults-result.json)。实际运行八条各 `0.25 km` 短轨迹：8/12 MHz 真空 O/X，以及平滑磁化 O/X 的规范 L2 与同进程显式冻结 L2 对照。真空位置最大误差约 `6.43e-13 km`，群程最大误差约 `1.39e-16 km`；磁化对照位置差为 0。替身包装器观测到 DOP853 消费的 `atol_k/k0=1e-12`。磁化对照检查参数路由及复现，不构成独立物理参考。
2. `python -B tasks/AH3D-PH474-L2-DEFAULT-20261103/bounded_checks.py physalign`：**PASS**，包括额外两条各 `0.05 km` 短轨迹；见[回执](physalign-receipt.json)和[结果](physalign/engineering_result.json)。该现行回归输出中保留了来源测试的历史 S4 任务标识；本次执行身份以当前目录与回执为准。
3. `main2.py --help` 退出 0，CLI 对旧 solve/final 档分别退出 2；见 [CLI 证据](cli-evidence.json)。路由测试覆盖 fan、P2P、默认与显式系数和同档复用；这不是完整 P2P 场景运行。
4. 完成上述完整短轨迹检查后，唯一源码补改是 `feature_gates` 中配置项缺省时的标签回退值，并在新回归中增加相应断言；未再改积分或参数计算。用 `runpy` 单独调用该回归的 `settings_checks()`，**PASS**，见[定点记录](metadata-settings.log)。因此八条轨迹的执行字节版本早于此标签补改，不宣称它们发生在补改之后。

## 边界与交接

未运行 PH、MATLAB、IRI、IGRF、REA、完整 P2P、完整生产场景或正式 26 例矩阵。短轨迹及路由结果不构成一般三维 O/X、自然端点、长轨迹性能或散射功率的科学验收。ECEF DOP853 核心静态检查未发现独立 accepted-step 数量帽；长轨迹成本仍待实际测量。

本次迁移会使绑定活动源码或文档的历史科学清单出现预期哈希差异；历史结果和旧快照不因此改变或自动迁移。总体/S5/S6 导航、历史清单版本映射及本包最终绑定交 Root 完成。本文为 Worker 自检和事实交接，**不是独立 ACCEPT**；最终结论由实际 Reviewer 审查。
