# AH 脚本与高精度路径索引

2026-11-03；HEAD `4be052c0c69f31a631beb8ed823c13d539036274`。本索引区分当前生产、精度选项和历史复核入口；不把脚本整理当作新增科学验收。

## 当前生产链

`main.py / main2.py → models/unified_hybrid_pipeline.py → ECEF AH3D / P2P → Appleton–Hartree dispersion`。

| 文件 | 保留用途 |
|---|---|
| [ray_tracer_ah_3d_ecef.py](./models/ray_tracer_ah_3d_ecef.py) | ECEF DOP853、完整 O/X Hamiltonian、真实磁场、投影和自然事件。 |
| [dispersion_ah.py](./models/dispersion_ah.py) | 色散分支、导数和群速度；连续性接口仍有其他调用，保留。 |
| [point_to_point_ah_3d_ecef.py](./models/point_to_point_ah_3d_ecef.py) | ECEF P2P 多尺度求解及最终高精度复核。 |
| [unified_hybrid_pipeline.py](./models/unified_hybrid_pipeline.py) | 正式请求、工厂和参数传递。 |
| [geomagnetic_field_3d.py](./models/geomagnetic_field_3d.py) | 真实 IGRF cubic，保留已修严格样条求解精度。 |
| [absorption_collision.py](./models/absorption_collision.py) | 消费射线实际等离子体系数的碰撞吸收。 |

## 精度与可用证据

[config.py](./config.py) 统一正式生产及研究入口为 `L2`：`rtol=1e-10`、`atol_r_km=1e-8`、`atol_p=1e-12`、步长 `0.002–0.5 km`，普通空间差分名义基准 `0.03125 km`。积分器按实际频率的 `k0` 将 `atol_p` 换算为物理 `atol_k`。fan、P2P 求解和最终结果都使用 L2；同档收敛轨迹可复用，结果记录是否回算。正式入口拒绝旧低档，`smoke`、`research`、`reference` 等仍供底层显式诊断。细结构可令差分更小；显式径向兼容模式保留其旧配置。

L2 会增加长轨迹和 P2P 的计算量。此次入口迁移只用有界短程验证参数及路由，不自动获得完整场景、一般三维、截止或 P2P 完整性资格。当前回归见 [`_verify_l2_defaults.py`](./tests/_verify_l2_defaults.py)，历史绑定脚本和数据不重跑或覆盖；旧清单绑定的活动源文件可能因已授权迁移发生预期哈希变化。

最新已接受的球面 IGRF O/X 反射细档为任务专属 L2：显式 SI、`spatial_fd_base_km=0.03125`，以及冻结的积分/事件/支持配置。其[采集入口和参数](./tasks/AH3D-PH474-S5-REFLECTIGRF-20261103/run.py)、[结果](./tasks/AH3D-PH474-S5-REFLECTIGRF-20261103-report.md)及[实际终审](./tasks/AH3D-PH474-S5-REFLECTIGRF-20261103-review.md)原位保留。该旧任务入口绑定其写入目录与预算，不应直接再次执行；新实验须使用新任务。WGS84、一般三维和截止退化资格不由此自动获得。

正式轻量检查入口是 [tests/_run_production_checks.py](./tests/_run_production_checks.py)，只执行它明确列出的七项；当前物理/SI回归与历史分类见[测试入口索引](./tests/README.md)。不要 glob 运行全部 `_verify*.py`。历史 AH 测试可能针对已隔离的 legacy 路径，并非当前生产验收。

## 历史归档

两份旧 V2 测试及其他历史脚本按[可回滚归档清单](./tasks/AH3D-PH474-AH-CLEANUP-20261103-archive-list.md)处理；归档方法选择和实际完成状态以该任务报告为准。冻结失败证据和必要负控制保留；历史清理当时未迁移生产默认，后续 L2 默认迁移另有授权，见[L2 迁移记录](./tasks/AH3D-PH474-L2-DEFAULT-20261103/worker-report.md)。
