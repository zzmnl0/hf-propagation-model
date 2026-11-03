# AH 历史脚本可回滚归档清单

状态：`LOGICAL_ARCHIVE_COMPLETED / ORIGINAL_PATHS_RETAINED`。用户选择“原位逻辑归档，保持冻结校验”；两份旧 V2 入口已从活跃使用索引中退出，原字节/路径不变。已保存恢复快照和 SHA，物理移动/删除0。依据为当前源码和实际 Worker 的只读调用/冻结引用盘点，不把静态不兼容描述为本轮运行 FAIL。

| 原路径 | 现状与用途 | 推荐动作 | 若选择物理归档的目标 |
|---|---|---|---|
| `tests/_verify_ah3d_v2_default.py` | 旧 legacy AH V2 默认验收，断言 solver=`v2`；当前正式 solver 为 `ecef_dop853`。被旧 S1 冻结清单按当前路径引用。 | 原位逻辑归档，退出活跃使用 | `legacy/tests/ah3d_v2_retired/_verify_ah3d_v2_default.py` |
| `tests/_verify_ah3d_v2_scenarios.py` | 旧 legacy AH V2 smoke 场景测试；同样被旧 S1 冻结。 | 原位逻辑归档，退出活跃使用 | `legacy/tests/ah3d_v2_retired/_verify_ah3d_v2_scenarios.py` |

精确原 SHA：

- `_verify_ah3d_v2_default.py`：`D73C09819763FFE724644E256F69F9DDE598F6352A6327C2493B5BF286739329`。
- `_verify_ah3d_v2_scenarios.py`：`FFB98F03AEC16748689F71985651B077FC234DCE2FBFADF2F240D51DDAC57A68`。

冻结证据位置：[S1 原清单](../output2/a4_pharlap/diag_s1_20261013_v1/manifest.json:2696)、[S1 续接清单](../output2/a4_pharlap/diag_s1_resume_20261014_v1/manifest.json:3428)。后者归类为 `review_evidence`。另存 `source_before/tests/` 副本，但副本不会自动替代清单中的原绝对路径。生产检查中的 `excluded_tests` 只是排除记录，不是运行记录。

按两个精确文件名，在 `tasks` 与 `output2/a4_pharlap` 的 JSON/Python 中包含 ignored 文件搜索，未发现 Python 按完整文件名调用；这不保证仓库外和动态调用不存在。

## 扩大盘点和其余文件处理

[逐测试冻结/调用索引](./AH3D-PH474-AH-CLEANUP-20261103/legacy-test-inventory.json)覆盖根目录152个Python及7个MATLAB文件。两份S1各直接绑定150/152个Python和7/7个MATLAB；仅不在旧S1中的两项为当前有效`_verify_ah3d_physalign.py`、`_verify_siaccess_request.py`，应保留。38项直接导入legacy的脚本全部绑定，其中一部分仍有共享检验/历史调用用途，不一律标“无用”。

历史绑定不等于当前字节全部相同：既有dirty文件`_verify_p3_compatibility.py`已在本轮前改变；新索引保留当前与历史SHA差异，不把本次归档说成修复了旧清单版本漂移。38项legacy直导入文件与两S1的SHA均一致。两份V2的归档仅改变索引角色，没有改源码或旧清单。

| 范围 | 动作及原因 |
|---|---|
| 正式 `models/` 中 ECEF AH3D、色散、P2P、统一入口、IGRF 与吸收 | 全部保留；当前高精度生产路径，没有证据支持删除关键符号。 |
| `config.py` 各积分档和必要粗档/负控制 | 保留；粗档用于性能和收敛诊断，不能按单次 FAIL 全删。生产默认不变。 |
| 121 份 AH 任务 Python 脚本 | [逐项原位索引](./AH3D-PH474-AH-CLEANUP-20261103/historical-script-index.json)，只作历史复核；本轮不移动。 |
| `output2/a4_pharlap` 的冻结脚本/raw/reference/manifest | 维持原路径与原字节；继续按各冻结清单解释，不删除失败证据。 |
| `legacy/models/` 三个 AH 模块 | 保留现有隔离；正式 source-separation 测试显式要求存在。 |
| `_verify_ah3d_dispersion.py`、`_verify_AH3D.py`、历史 PH 冻结入口 | 保留；有共享检验或调用/冻结引用，不能按名称和历史 FAIL 判断冗余。 |

## 回滚及物理归档影响

两份候选原字节已存于[恢复映射](./AH3D-PH474-AH-CLEANUP-20261103/rollback-manifest.json)所列 `baseline/tests/`。若选物理归档，归档前后核 SHA 和长度；恢复时先确认原路径无其他新修改，再从核过 SHA 的快照恢复。归档副本只是历史源码，不承诺移位后直接执行；复核旧结果须恢复原路径及对应历史配置。

物理移动会使依赖原绝对路径的旧 manifest 校验在恢复前失败。推荐原位逻辑归档可以保持这一证据链连续。两种方式均不涉及生产数值代码、有效高精度链或科学门变化。
