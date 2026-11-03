# 当前 AH 测试入口与历史归档

正式生产与高精度研究路径见[AH脚本索引](../AH3D-SCRIPTS.md)。本索引不运行科学，也不把框架检查当物理验收。

## 正式轻量检查

[`_run_production_checks.py`](./_run_production_checks.py)明确选择七项：

- `_verify_production_boundary.py`
- `_verify_a4_c0r1_geometry.py`
- `_verify_a4_c0r3_ecef_hamiltonian.py`
- `_verify_absorption_consistency.py`
- `_verify_a4_c0r6_entry_migration.py`
- `_verify_legacy_source_separation.py`
- [`_verify_l2_defaults.py`](./_verify_l2_defaults.py)：L2 参数、正式路由和短程轨迹回归

该完整集合含受控 IRI 原生调用，不应把列出测试当作已执行；可用 `--test _verify_l2_defaults.py` 只运行当前 L2 检查。额外当前回归为[物理系数及吸收一致性](./_verify_ah3d_physalign.py)。高精度轨迹的 FD、积分、支持和事件资格使用独立任务契约，不由文件名或轻量测试通过推定。

## 已退出活跃使用的旧入口

`_verify_ah3d_v2_default.py`、`_verify_ah3d_v2_scenarios.py`属于旧legacy V2路线，原位逻辑归档，不作为现行ECEF/DOP853研究或生产验收命令。`_verify_siaccess_request.py` 绑定历史任务的旧档路由和调用契约，原样保留供旧版本复核；当前统一 L2 的请求透传由 `_verify_l2_defaults.py` 检查。

两份S1清单各绑定根目录150/152个Python和全部7个MATLAB历史入口，因此保留这些原路径/字节。38项直接legacy导入不等于全无用；其中历史色散、Jacobian、参考和P5复核仍有共享/历史调用。需要复核时按各自冻结配置运行，禁止glob全部测试或将旧FAIL覆盖为当前PASS。

逐项用途及历史冻结值见[机器清单](../tasks/AH3D-PH474-AH-CLEANUP-20261103/legacy-test-inventory.json)，两份旧V2的快照与恢复方式见[归档清单](../tasks/AH3D-PH474-AH-CLEANUP-20261103-archive-list.md)。该清理任务当时未做默认迁移；后续已授权的 L2 迁移按新任务记录和 before 快照追溯。
