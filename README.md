# HF 短波混合传播模型 — 详细使用手册

**研究项目**：李慧敏-22  复杂电离层环境短波反射/散射混合传播机理及建模研究

**测试链路**：TX @ 30°N 120°E → RX @ 1169 km，f = 10 MHz，2020-06-01 12:00 UT

**已完成阶段**：Phase 1–7 全部完成，全部 9 个场景测试通过，无卡顿

**当前公式修复状态（2026-07-01）**：已修正群路径/群延迟为 `integral(ds/n)`、Appleton-Hartree O/X 分支、3D D 层吸收、3D Born/Rytov/MPS 诊断耦合，以及 P2P_3D forward/central Jacobian 选项。

---

## 目录

1. [快速上手](#1-快速上手)
2. [目录结构](#2-目录结构)
3. [整体数据流](#3-整体数据流)
4. [模块详解](#4-模块详解)
5. [参数速查表](#5-参数速查表)
6. [输出格式](#6-输出格式)
7. [场景运行示例](#7-场景运行示例)
8. [验证脚本汇总](#8-验证脚本汇总)
9. [常见修改方法](#9-常见修改方法)
10. [性能参考](#10-性能参考)
11. [参考文献](#11-参考文献)

---

## 1. 快速上手

```bash
# 环境（仅支持 iri2016，不能用 iricore）
conda activate pytorch_cpu

# 完整运行所有场景 -> output/
conda run -n pytorch_cpu python main.py

# 单独场景
conda run -n pytorch_cpu python -c "from main import run_baseline; run_baseline()"
conda run -n pytorch_cpu python -c "from main import run_tube_radar; run_tube_radar()"

# 背景电子密度图
conda run -n pytorch_cpu python viz/plot_ne_background.py

# 验证各阶段（全部应输出 0 fail(s)）
conda run -n pytorch_cpu python tests/_verify_part1.py
conda run -n pytorch_cpu python tests/_verify_part6.py
conda run -n pytorch_cpu python tests/_verify_group_path.py
conda run -n pytorch_cpu python tests/_verify_absorption_3d.py
conda run -n pytorch_cpu python tests/_verify_scatter_coupling.py
conda run -n pytorch_cpu python tests/_verify_p2p_jacobian.py
conda run -n pytorch_cpu python tests/_verify_tube.py
conda run -n pytorch_cpu python tests/_verify_phase5.py
conda run -n pytorch_cpu python tests/_verify_phase6.py
conda run -n pytorch_cpu python tests/_verify_phase7.py
```

**依赖**：`numpy scipy matplotlib iri2016 ppigrf`

> **Windows 编码**：`print()` 中只用 ASCII；中文注释不影响运行。

---

## 2. 目录结构

```
code/
├── config.py                     # 全局参数（唯一配置入口）
├── utils.py                      # 物理工具函数（EM、功率、坐标、地磁）
├── main.py                       # 场景驱动脚本（9 个场景函数）
│
├── models/
│   ├── ionosphere_model.py       # M1:  IRI + TID / Es / 等离子体泡 / Spread-F
│   ├── geomag_field.py           # B2:  三维 IGRF（ppigrf，Phase 5）
│   ├── spread_f_model.py         # M3b: Rino 1979 幂律相位屏
│   ├── ray_tracer.py             # M2a: 二维 Haselgrove RK4 + RefractiveIndexAH
│   ├── point_to_point.py         # M2b: 二维 P2P 变分求解（Nosikov 2020）
│   ├── tube_tracer.py            # A3:  二维 flux-tube（Coleman 1997/1998）
│   ├── ray_tracer_3d.py          # M2d: 三维 Haselgrove + TID-3D + 场向泡（Phase 6/7）
│   ├── point_to_point_3d.py      # M2e: 三维 P2P Newton 双参数 shooting（Phase 6）
│   ├── tube_tracer_3d.py         # A3b: 三维 flux-tube 椭圆截面（Phase 7）
│   ├── es_model.py               # M3:  Hao 2017 三段式 + Booker-Gordon（Phase 7）
│   ├── pe_propagator.py          # M4:  PE/SSF 等离子体泡散射（Carrano 2020）
│   ├── scatter_selector.py       # C0:  Born/Rytov/MPS 方法选择器（Phase 6）
│   ├── scatter_born.py           # C1:  Born 弱散射体散射截面（Phase 6）
│   ├── scatter_rytov.py          # C2:  Rytov 闪烁理论（Phase 6）
│   ├── scatter_coupler.py        # C6:  三维模式 Born/Rytov/MPS 统一诊断耦合
│   ├── scatter_phase_screen.py   # C3/C5: Rino 二维相位屏 + Knepp MPS（Phase 7）
│   └── hybrid_model.py           # M5-M7: 完整 pipeline（P-D 谱 + 主模式）
│
├── viz/
│   ├── plot_utils.py             # 共享绘图函数
│   └── plot_ne_background.py     # 4 面板 Ne 背景图
│
├── analysis/
│   ├── compare_pd.py             # 实测 vs 模型 P-D 谱比对（CLI）
│   └── mode_analyzer.py          # 模式特征分析工具（CLI）
│
├── tests/
│   ├── _verify_part0.py          # 导入健全性
│   ├── _verify_part1.py          # IonosphereModel 验证
│   ├── _verify_part2.py          # 二维射线追踪验证
│   ├── _verify_part3.py          # P2P 变分验证
│   ├── _verify_part4.py          # Es 模型验证
│   ├── _verify_part5.py          # PE/SSF 验证
│   ├── _verify_part6.py          # 混合 pipeline 验证
│   ├── _verify_radar.py          # OTH 雷达几何验证
│   ├── _verify_OX.py             # O/X 磁离子分裂验证
│   ├── _verify_absorption.py     # D 层吸收验证
│   ├── _verify_absorption_3d.py  # 三维功率链 D 层吸收验证
│   ├── _verify_group_path.py     # 几何路径/群路径/群延迟公式验证
│   ├── _verify_p2p_jacobian.py   # 三维 P2P forward/central Jacobian 验证
│   ├── _verify_scatter_coupling.py # 三维 Born/Rytov/MPS 耦合验证
│   ├── _verify_spreadf.py        # Spread-F 相位屏验证
│   ├── _verify_tube.py           # 二维流管追踪验证
│   ├── _verify_phase5.py         # Phase 5：三维背景 + 多方向 TID 验证
│   ├── _verify_phase6.py         # Phase 6：三维射线 + Born/Rytov 验证
│   └── _verify_phase7.py         # Phase 7：泡/流管/BG/屏/MPS 验证
│
├── upgrade_plan_3d_scatter.txt   # Phase 5/6/7 技术路线图
└── output/                       # 生成图 / CSV（gitignored）
```

---

## 3. 整体数据流

### 3.1 二维路径（主线，Phase 1–4）

```
config.py
    |
    v
IonosphereModel.build_Ne_field(BG_X, BG_Z)   [M1]
    ├─ IRI-2016 背景 profile  (iri2016.IRI)
    ├─ TID 扰动               (Hooke 1968, enable=TID['enable'])
    ├─ Es 薄层                (Hao 2017,   enable=ES['enable'])
    ├─ 等离子体泡             (Gaussian,   enable=BUBBLE['enable'])
    └─ Spread-F               (Rino 1979,  enable=SPREAD_F['enable'])
    |
    v  Ne_2d (Nx, Nz) [m^-3]
    |
    v
RefractiveIndex(Ne_2d, BG_X, BG_Z, FREQ_MHZ)           各向同性  [M2a]
RefractiveIndexAH(Ne_2d, BG_X, BG_Z, FREQ_MHZ, 'O'/'X') Appleton-Hartree
    |
    ├──[tube_mode=False]──────────────────────────────────────────────────
    |   find_all_rays_p2p(tx, rx, n_model, freq) -> rays [M2b]
    |       光程泛函 S = sum n(mid)*|seg|，梯度下降 + 弹簧平滑
    |
    └──[tube_mode=True]───────────────────────────────────────────────────
        TubeRayTracer.compute(tx, x_tgt, Pt, Gt, Gr, sigma0) -> modes [A3]
            扇形扫射 -> 相邻射线对 -> Newton 精化 -> Coleman 功率公式
    |
    v
Power pipeline（每条 ray）:
    Friis / 雷达方程 / Coleman 流管  ->  D层修正  ->  Es修正  ->  PE/SSF修正
    |
    v
build_pd_spectrum() -> tau_ax, pd_W
identify_main_mode() -> main_mode, ranked_modes, summary
```

### 3.2 三维路径（Phase 5–7）

```
IonosphereModel.build_Ne_field()  (同二维，建 Ne_2d 背景)
    |
    v
RefractiveIndex3D(Ne_2d, BG_X, BG_Z, FREQ_MHZ,
                  tid_params=TID,            # 三维 TID（ky != 0）
                  bubble_3d_params=BUBBLE_3D, # 场向椭球耗空
                  earth_flat=True)           # Jones 1975 展平

    n2(x, y, z):
        z_phys = z_cart - (x^2+y^2)/(2*RE)
        Ne_tot = Ne_bg + dNe_TID - Ne_bubble  （可叠加）
        n2     = max(1 - fp^2/f^2, 1e-6) * (1+z_phys/RE)^2
    |
    ├──[单射线]──────────────────────────────────────────────────────────
    |   trace_single_ray_3d(tx, beta_deg, az_deg, n3d) -> ray_3d [M2d]
    |       state = [x,y,z,kx,ky,kz], Haselgrove RK4, on-shell renorm
    |
    ├──[P2P]────────────────────────────────────────────────────────────
    |   find_all_rays_3d(tx, rx, n3d, freq, p2p_params=P2P_3D) [M2e]
    |       Newton(beta,az)，2x2 Jacobian 可选 forward / central 差分
    |
    └──[流管]───────────────────────────────────────────────────────────
        TubeRayTracer3D.compute(tx, x_tgt, Pt, Gt, Gr, sigma0) [A3b]
            F_focus_3d = sqrt(F_focus_el * F_focus_az)
            A_tube_3d  = delta_x * delta_y
    |
    v
散射模块（Phase 6/7，按 scatter_selector 选择）:
    scatter_phase_variance() -> sigma_phi^2
        < 0.1   -> Born  (scatter_born.py)
        0.1-10  -> Rytov (scatter_rytov.py)
        >= 10   -> MPS   (scatter_phase_screen.py)
    scatter_coupler.apply_scatter_regime_3d()
        给三维 mode dict 添加 scatter_method / sigma_phi / Born 或 Rytov 诊断字段

Es 各向异性散射（Phase 7）:
    booker_gordon_sigma_v(freq, Cs, Ne, L_perp, L_para, p_es, aspect)
        von Karman 谱 (perp) x Gaussian 方位角滤波器 (para)
```

---

## 4. 模块详解

### M1 — IonosphereModel

**文件**：`models/ionosphere_model.py`

```python
from models.ionosphere_model import IonosphereModel
import config as cfg

iono = IonosphereModel(
    iri_params      = {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
    tid_params      = {**cfg.TID,      'enable': True},   # 可选
    es_params       = {**cfg.ES,       'enable': True},   # 可选
    bubble_params   = {**cfg.BUBBLE,   'enable': True},   # 可选
    spread_f_params = {**cfg.SPREAD_F, 'enable': True},   # 可选
    freq_MHz        = cfg.FREQ_MHZ,
)
Ne_2d, n_2d = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z, t=0.0)
# Ne_2d: (Nx, Nz) [m^-3]   n_2d: (Nx, Nz) 各向同性折射率
```

**叠加顺序**（各层独立开关，顺序固定）：

| 顺序 | 扰动 | 物理模型 | 配置组 |
|------|------|----------|--------|
| 1 | IRI 背景 | iri2016，路径中点单剖面或每 50 km 横向采样 | `IRI_DT/LAT/LON`, `IRI_LATERAL` |
| 2 | 多方向 TID | Hooke (1968)，N 个方向分量叠加 | `TID` |
| 3 | Es 薄层 | Hao (2017) 密度剖面 | `ES` |
| 4 | 等离子体泡（二维） | Gaussian 耗空 (x,z) 平面 | `BUBBLE` |
| 5 | 等离子体泡（三维） | 场向椭球（RefractiveIndex3D 中处理） | `BUBBLE_3D` |
| 6 | Spread-F | Rino (1979) 幂律一维随机场 | `SPREAD_F` |

---

### M2a — 二维 RefractiveIndex / RefractiveIndexAH

**文件**：`models/ray_tracer.py`

```python
from models.ray_tracer import RefractiveIndex, RefractiveIndexAH
from models.ray_tracer import trace_single_ray, shoot_rays_fan

# 各向同性
nm = RefractiveIndex(Ne_2d, BG_X, BG_Z, FREQ_MHZ)

# O/X 磁离子
nm_O = RefractiveIndexAH(Ne_2d, BG_X, BG_Z, FREQ_MHZ,
                          wave_mode='O', geomag=cfg.GEOMAG)
nm_X = RefractiveIndexAH(Ne_2d, BG_X, BG_Z, FREQ_MHZ, wave_mode='X', geomag=cfg.GEOMAG)

# 单条射线
ray = trace_single_ray(TX_POS, beta_deg=25.0, n_model=nm, freq_MHz=FREQ_MHZ)
# ray 字典键: beta_deg, trajectory [(x,z),...], geom_path_km, group_path_km, tau_ms,
#             h_reflect_km, beta_recv_deg, wave_mode, label

# 扇形扫射
rays = shoot_rays_fan(TX_POS, nm, rt_params=cfg.RT)
```

**Haselgrove 二维方程**（等效于 Jones & Stephenson 1975）：

```
dx/ds = px/|p|,   dz/ds = pz/|p|
dpx/ds = (k0^2/2|p|) * dn^2/dx
dpz/ds = (k0^2/2|p|) * dn^2/dz
|p| = n * k0 (on-shell 约束)
```

**路径/时延约定**：

- `geom_path_km = integral(ds)`：射线几何弧长。
- `group_path_km = integral(ds/n)`：冷等离子体近似下的群路径（`n_g ~= 1/n`）。
- `tau_ms = group_path_km / c`，所有功率链中的群路径字段均按此口径使用。

---

### M2b — 二维 P2P 变分求解（Nosikov 2020）

**文件**：`models/point_to_point.py`

```python
from models.point_to_point import find_all_rays_p2p

rays = find_all_rays_p2p(
    tx_km      = TX_POS,
    rx_km      = RX_POS,
    n_model    = nm,
    freq_MHz   = FREQ_MHZ,
    p2p_params = cfg.P2P,
    wave_mode  = None,      # None=各向同性 | 'O' | 'X' | 'both'
    geomag     = cfg.GEOMAG,
)
```

**模式标签分类规则**：

| h_reflect | 标签 | 变体（O/X 开启后） |
|-----------|------|------------------|
| < 140 km | `Es` | `Es_O` / `Es_X` |
| 140–200 km | `E` | `E_O` / `E_X` |
| 200–300 km, tau<5 ms | `1F_low` | `1F_low_O` / ... |
| 200–300 km, tau>=5 ms | `1F_high` | `1F_high_O` / ... |
| >= 300 km | `2F` | `2F_O` / `2F_X` |

---

### M2c — 二维 TubeRayTracer（Coleman 1997/1998）

**文件**：`models/tube_tracer.py`

```python
from models.tube_tracer import TubeRayTracer
import numpy as np

sigma0 = 10.0 ** (cfg.RADAR['sigma0_ground_dB'] / 10.0)
tracer = TubeRayTracer(n_model=nm, freq_MHz=FREQ_MHZ,
                        tube_params=cfg.TUBE_TRACER, rt_params=cfg.RT)
modes = tracer.compute(TX_POS, x_tgt=1169.0,
                       Pt_W=1000.0, Gt=1.0, Gr=1.0, sigma0=sigma0, newton=True)
# 额外字段: F_focus, A_tube_km2, tau_spread_ms, tau_2way_ms
```

**核心公式（Coleman 1997, Eq. 15-17）**：

```
F_focus   = R_eff * d_beta_rad / d_x_land    # 聚焦因子（1=自由空间）
A_tube    = delta_x * L_cross_km              # 地面投影面积 [km^2]
Pr        = Pt*Gt*Gr*lambda^2*sigma0*A_tube / ((4*pi)^3 * P_one^4)
f_pulse   = T_pulse / tau_spread  if tau_spread > T_pulse else 1.0
```

---

### M2d — 三维 RefractiveIndex3D 与射线追踪（Phase 6/7）

**文件**：`models/ray_tracer_3d.py`

```python
from models.ray_tracer_3d import RefractiveIndex3D, trace_single_ray_3d
from models.ray_tracer_3d import shoot_rays_fan_3d

# 创建三维折射率模型
n3d = RefractiveIndex3D(
    Ne_2d            = Ne,
    x_km             = BG_X,
    z_km             = BG_Z,
    freq_MHz         = FREQ_MHZ,
    tid_params       = {**cfg.TID,       'enable': True},   # 含 ky 分量的三维 TID
    bubble_3d_params = {**cfg.BUBBLE_3D, 'enable': True},   # 场向椭球耗空
    earth_flat       = True,   # Jones 1975 展平
)

# 单条三维射线
ray3d = trace_single_ray_3d(TX_POS, beta_deg=25.0, az_deg=5.0,
                              n3d=n3d, freq_MHz=FREQ_MHZ, rt_params=cfg.RT_3D)
# 额外字段: x_land_km, y_land_km, z_land_km, h_reflect_km,
#           beta_deg_tx, az_deg_tx, beta_deg_rx, az_deg_rx,
#           azimuth_deflect_deg, trajectory_3d [(x,y,z,kx,ky,kz),...]
```

**地球展平变换（Jones 1975）**：

```
z_phys = z_cart - (x^2+y^2)/(2*RE)
n2_eff(x,y,z_cart) = n2_plasma(x, z_phys) * (1 + z_phys/RE)^2
```

**三维 TID 扰动（Phase 5，Hooke 1968 多方向）**：

```
dNe(x,y,z,t) = sum_i A_i * [-dNe_dz*sinI*sin(phi_i) - k_para_i*Ne*cos(phi_i)]
phi_i = kx_i*x + ky_i*y + kz_i*z - omega_i*t
kx_i = k_h*cos(az_i - link_az),  ky_i = k_h*sin(az_i - link_az)
```

**三维场向等离子体泡（Phase 7，Ma et al. 2023）**：

```
r_perp = r - (r . B_hat) * B_hat     # 垂直磁场方向
r_para = r . B_hat                    # 沿磁场方向
frac   = delta_max * exp(-|r_perp|^2/sigma_perp^2 - r_para^2/sigma_para^2)
Ne_tot = Ne_bg * (1 - frac)
B_hat  = (cos(I)*cos(D-phi), cos(I)*sin(D-phi), -sin(I))   在(x,y,z)系
```

---

### M2e — 三维 P2P Newton 双参数 Shooting（Phase 6）

**文件**：`models/point_to_point_3d.py`

```python
from models.point_to_point_3d import find_all_rays_3d

modes_3d = find_all_rays_3d(
    tx_pos     = TX_POS,
    rx_pos     = RX_POS,
    n3d        = n3d,
    freq_MHz   = FREQ_MHZ,
    p2p_params = cfg.P2P_3D,
)
# 额外字段（比二维 P2P 新增）:
#   y_land_km, azimuth_deflect_deg, beta_deg_tx, az_deg_tx
```

**Newton 迭代（2×2 Jacobian，可选 forward / central 差分）**：

```
F1(beta, az) = x_land - x_tgt = 0
F2(beta, az) = y_land - y_tgt = 0
J = [[dF1/dbeta, dF1/daz],      # 2x2 Jacobian，步长 d_beta=0.05°, d_az=0.10°
     [dF2/dbeta, dF2/daz]]
[beta, az] -= J^{-1} * F
```

`P2P_3D['jacobian']` 默认 `'forward'`（每轮 3 条射线，速度快），可改为 `'central'`（每轮 5 条射线，精度更高）。

---

### A3b — 三维 TubeRayTracer3D（Phase 7）

**文件**：`models/tube_tracer_3d.py`

```python
from models.tube_tracer_3d import TubeRayTracer3D

tt3d = TubeRayTracer3D(n3d=n3d, freq_MHz=FREQ_MHZ,
                        tube_params=cfg.TUBE_3D, rt_params=cfg.RT_3D)
modes3d = tt3d.compute(TX_POS, x_tgt=1169.0,
                        Pt_W=1000.0, Gt=1.0, Gr=1.0, sigma0=sigma0)
# 额外字段: F_focus_el, F_focus_az, F_focus_3d, A_tube_km2
```

**三维聚焦公式（Coleman 1997，扩展到三维）**：

```
F_focus_el = R_eff * d_beta_rad / |d_x_land|          # 仰角平面聚焦
F_focus_az = R_eff * d_az_rad  / |d_y_land(+az)-d_y_land(-az)|  # 方位角平面聚焦
F_focus_3d = sqrt(F_focus_el * F_focus_az)             # 三维几何平均
A_tube_3d  = |d_x_land| * |d_y_land|                  # 椭圆截面面积 [km^2]
```

---

### M3 — Es 层模型（Hao 2017 + Booker-Gordon，Phase 7）

**文件**：`models/es_model.py`

**① Hao (2017) 三段式模型（类 EsLayerModel）**：

```python
from models.es_model import EsLayerModel

es = EsLayerModel(foEs_MHz=5.0, h_Es_km=110.0, delta_h_m=115.0, n_exp=5,
                   L1_m=300.0, L2_m=300.0, L3_m=30.0, delta_N_N=0.3,
                   fr=0.25, fs=0.10)
result = es.compute_power(Pt_W=1e3, Gt=1.0, Gr=1.0, f_MHz=10.0,
                           D_km=1200.0, theta_rad=np.radians(20.0))
# result 键: Pr_W, Pr_reflect_W, Pr_scatter_W, mode, alpha, tau_ms, delta_tau_ms
```

| foEs/f 区间 | 模式 | 反射权重 alpha |
|-------------|------|---------------|
| > fr=0.25 | reflect | 1.0 |
| fs..fr | mixed | 线性插值 |
| < fs=0.10 | scatter | 0.0 |

**② Booker-Gordon 各向异性散射（Phase 7，函数级 API）**：

```python
from models.es_model import (booker_gordon_sigma_v,
                               es_aspect_angle_deg,
                               es_scatter_valid)

# 最优方位角（垂直 B，aspect=0）时的体散射截面
sv = booker_gordon_sigma_v(
    freq_MHz   = 10.0,
    Cs_rel     = 0.10,           # 相对起伏幅度
    Ne_m3      = 2e10,           # 背景电子密度 [m^-3]
    L_perp_m   = cfg.ES_BG['L_perp_m'],   # 外尺度（垂直 B）[m]
    L_para_m   = cfg.ES_BG['L_para_m'],   # 场向展伸尺度 [m]
    p_es       = cfg.ES_BG['p_es'],        # 谱指数
    aspect_deg = 0.0,            # 方位角（0=最优散射）
)

# 根据几何计算有效方位角
asp = es_aspect_angle_deg(theta_rad=np.radians(80.0))   # 近垂直入射

# 检查方位条件是否满足
ok = es_scatter_valid(aspect_deg=asp, L_para_m=1e4, freq_MHz=10.0)
```

**Booker-Gordon 公式**（von Karman 谱 + Gaussian 方位角滤波器）：

```
sigma_v = r_e^2 * (2pi)^3 * C_s_abs * (Dk_perp^2 + kappa_0^2)^(-(p+3)/2)
          * exp(-Dk_para^2 * L_para^2 / 4)

C_s_abs = Cs_rel^2 * Ne^2 * p * (2pi)^2 * kappa_0^p    [Born 惯例]
kappa_0 = 2*pi / L_perp      Dk = 2*k0（反散射 Bragg 波数）
Dk_perp = Dk*cos(asp),  Dk_para = Dk*sin(asp)
```

---

### C0/C1/C2 — 散射选择器 + Born + Rytov（Phase 6）

**文件**：`scatter_selector.py`, `scatter_born.py`, `scatter_rytov.py`

```python
from models.scatter_selector import scatter_phase_variance, select_scatter_method
from models.scatter_born import born_sigma_v, born_scatter
from models.scatter_rytov import rytov_full

# 步骤 1：计算相位方差选择方法
sp = scatter_phase_variance(Ne0_m3=5e10, Cs_rel=0.01,
                             DeltaL_km=100.0, L_outer_km=10.0,
                             freq_MHz=10.0)   # -> sigma_phi^2 [rad^2]
method = select_scatter_method(sp)            # -> 'Born'/'Rytov'/'MPS'

# Born 弱散射（sigma_phi^2 < 0.1）
sv = born_sigma_v(freq_MHz=10.0, Cs_rel=0.01, p_spec=3.0,
                   L_outer_km=10.0, Ne0_m3=5e10)   # sigma_v [m^-1]

# Rytov 中等散射（0.1 <= sigma_phi^2 < 10）
ryt = rytov_full(freq_MHz=10.0, Cs_rel=0.01, p=3.0,
                  L_outer_km=10.0, Ne0=5e10,
                  z_eff_km=300.0, DeltaL_km=100.0)
# ryt 键: S4, sigma_phi_rad, sigma_tau_ms, fresnel_r_km
```

**方法阈值（Rino 1979 / Yeh & Liu 1982）**：

| sigma_phi^2 [rad^2] | 方法 | 适用条件 |
|---------------------|------|---------|
| < 0.1 | Born | 弱单散射 dNe/Ne < 5% |
| 0.1 – 10 | Rytov | 中等散射 S4 < 0.6 |
| >= 10 | MPS | 强散射（等离子体泡） |

---

### C3/C5 — Rino 二维相位屏 + Knepp MPS（Phase 7）

**文件**：`models/scatter_phase_screen.py`

```python
from models.scatter_phase_screen import RinoPhaseScreen2D, MultiPhaseScreen

# ── C3: Rino (1979) 二维幂律相位屏 ──────────────────────────────────────────
screen = RinoPhaseScreen2D(
    freq_MHz   = 10.0,
    Cs_rel     = cfg.SCATTER['Cs_rel'],
    p          = cfg.SCATTER['p'],
    L_outer_km = cfg.SCATTER['L_outer_km'],
    Ne_m3      = 1e11,
    DeltaL_km  = cfg.SCATTER['DeltaL_km'],
    n_pts      = cfg.PHASE_SCREEN['n_pts'],    # 512
    dx_m       = cfg.PHASE_SCREEN['dx_m'],     # 200 m
    seed       = 42,
)
phi_screen = screen.generate()        # (n_pts, n_pts) [rad]
sigma_rms  = screen.sigma_phi_rms()   # 理论 RMS 相位 [rad]

# ── C5: Knepp (1983) 多相屏传播 ──────────────────────────────────────────────
mps = MultiPhaseScreen(
    freq_MHz   = 10.0,
    Cs_rel     = 0.10,
    p          = 3.0,
    L_outer_km = 10.0,
    Ne_m3      = 1e11,
    DeltaL_km  = 100.0,
    N_screen   = cfg.PHASE_SCREEN['N_screen'],  # 10
    n_pts      = 64,       # 快速测试用 64，精确用 512
    dx_m       = 200.0,
    seed       = None,     # None=随机
)
u0  = np.ones((64, 64), dtype=complex)   # 平面波入射
res = mps.propagate(u0)
# res 键: u_out (复场), I_mean, S4, Pr_norm（能量守恒比值）
```

**Rino 二维 PSD**（Rino 1979, Eq. 14）：

```
Phi_phi(kx,ky) = (r_e*lambda)^2 * C_s_abs * DeltaL * F(p)
                 * (kx^2 + ky^2 + kappa_0^2)^(-(p+2)/2)
F(p) = sqrt(pi) * Gamma((p+2)/2) / Gamma((p+3)/2)
sigma_phi^2 = A * kappa_0^(-p) / (2*pi*p)
```

**Knepp MPS 每步（Knepp 1983, Eq. 2-4）**：

```
Step 1: u -> u * exp(i*phi_i(x,y))          相位屏
Step 2: U = FFT2(u)                          角谱
Step 3: U -> U * exp(-i*(kx^2+ky^2)*dz/(2*k0))   自由传播
Step 4: u = IFFT2(U)
```

---

### M4 — PE/SSF 等离子体泡散射（Carrano 2020）

**文件**：`models/pe_propagator.py`

```python
from models.pe_propagator import PEPropagator, construct_incident_field

pe = PEPropagator(freq_MHz=FREQ_MHZ, pe_params=cfg.PE)
Ne_pe, x_pe, z_pe = pe.extract_domain(Ne_2d, BG_X, BG_Z,
                                        x_range=[500, 700], z_range=[200, 500])
u_inc = construct_incident_field(nm, ray, FREQ_MHZ, w0_km=20.0, x_pe, z_pe)
out   = pe.propagate(Ne_pe, x_pe, z_pe, u_inc, x_out=700.0, z_range_out=[200,500])
```

**SSF 单步**：

```
折射步（空域）: u_half = u * exp(i*k0*(n-1)*dx)
衍射步（谱域）: FFT -> exp(-i*kx^2*dx/(2k0)) -> IFFT
```

---

### M5-M7 — HybridPropagationModel（顶层 pipeline）

**文件**：`models/hybrid_model.py`

```python
from models.hybrid_model import HybridPropagationModel, build_pd_spectrum

model = HybridPropagationModel(
    iono_params  = {
        'iri_params':      {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
        'tid_params':      {**cfg.TID,      'enable': True},
        'es_params':       {**cfg.ES,       'enable': True},
        'bubble_params':   {**cfg.BUBBLE,   'enable': True},
        'spread_f_params': {**cfg.SPREAD_F, 'enable': True},
    },
    radar_params = {'freq_MHz': 10.0, 'Pt_W': 1e3,
                    'Gt': 1.0, 'Gr': 1.0,
                    'sigma0_ground_dB': -20.0},  # 流管模式必须
    geomag_params      = {**cfg.GEOMAG, 'enable_OX': True},   # O/X 分裂
    absorption_params  = {**cfg.ABSORPTION, 'enable': True},   # D 层修正
    radar_mode = True,    # 雷达功率公式
    tube_mode  = True,    # Coleman 流管（需 radar_mode=True）
)
modes, tau_ax, pd_W, main = model.compute(cfg.TX_POS, cfg.RX_POS, t=0.0)
summary = model.mode_summary   # 见 §6
```

---

## 5. 参数速查表

所有参数集中于 `config.py`，修改该文件即可全局生效。

### 5.1 测试链路基本参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `TX_LAT/TX_LON` | 30.0 / 120.0 | 发射机纬经度 [°] |
| `RX_RANGE` | 1169.0 | TX-目标水平距离 [km] |
| `FREQ_MHZ` | 10.0 | 工作频率 [MHz] |
| `PT_W` | 1000.0 | 发射功率 [W] |
| `GT`, `GR` | 1.0 | 天线增益（线性，各向同性） |
| `IRI_DT` | 2020-06-01 12:00 | IRI 快照时间 |
| `IRI_LAT/IRI_LON` | 32.5 / 120.0 | IRI 路径中点 [°N / °E] |
| `LINK_BEARING_DEG` | 0.0 | 链路方位角 [° 从正北顺时针] |

### 5.2 背景网格（`BG_X`, `BG_Z`）

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `BG_X_MIN/MAX/DX` | -50/1300/5 km | 水平范围和步长 |
| `BG_Z_MIN/MAX/DZ` | 60/600/2 km | 高度范围和步长 |

> 修改链路距离时同步更新 `BG_X_MAX = RX_RANGE + 200`

### 5.3 TID 参数（`config.TID`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | 总开关 |
| `amplitude` | 0.10 | 单分量峰值 dNe/Ne₀ |
| `lambda_h_km` | 300.0 | 水平波长 [km] |
| `T_s` | 2400.0 | 周期 [s]（40 min MSTID） |
| `I_dip_deg` | 50.0 | 地磁倾角 [°] |
| `H_km` | 60.0 | Chapman 标高 [km] |
| `omega_b_rad_s` | 2π/1200 | Brunt-Vaisala 频率 [rad/s] |
| `n_components` | 1 | 叠加分量数（Phase 5 多方向） |
| `az_deg_list` | [0.0] | 各分量传播方位角列表 [°] |
| `amplitude_list` | [0.10] | 各分量幅度列表 |
| `period_s_list` | [2400.0] | 各分量周期列表 [s] |
| `lambda_h_km_list` | [300.0] | 各分量水平波长列表 [km] |
| `link_bearing_deg` | 0.0 | 链路方位角（TID 投影用） |

### 5.4 Es 层参数（`config.ES`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | 总开关 |
| `foEs_MHz` | 5.0 | Es 等离子频率 [MHz] |
| `h_Es_km` | 110.0 | Es 中心高度 [km] |
| `delta_h_m` | 115.0 | 半厚度 [m]（Hao 2017 最佳拟合） |
| `n_exp` | 5 | 密度剖面指数 |
| `L1_m, L2_m` | 300.0 | 水平不规则体尺度 [m] |
| `L3_m` | 30.0 | 垂直不规则体尺度 [m] |
| `delta_N_N` | 0.3 | 相对密度起伏 ΔN/N |
| `fr` | 0.25 | 反射阈值 foEs/f |
| `fs` | 0.10 | 散射阈值 foEs/f |

### 5.5 Booker-Gordon Es 参数（`config.ES_BG`，Phase 7）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `L_perp_m` | 300.0 | 垂直 B 方向外尺度 [m]（von Karman 截断） |
| `L_para_m` | 10000.0 | 沿 B 方向尺度 [m]（方位角滤波器半宽） |
| `p_es` | 1.5 | Es 谱指数（浅于 F 层 p~3） |

> 方位角允许窗口：asp < arcsin(1/(Dk*L_para)) ≈ 0.014° @10MHz

### 5.6 等离子体泡参数

**二维 Gaussian（`config.BUBBLE`）**：

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | 开关 |
| `delta_max` | 0.6 | 最大耗散比（0–1） |
| `x0_km, z0_km` | 600.0 / 350.0 | 泡中心 [km] |
| `Lx_km, Lz_km` | 100.0 / 150.0 | 水平/垂直半宽 [km] |

**三维场向椭球（`config.BUBBLE_3D`，Phase 7）**：

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | 开关（在 RefractiveIndex3D 中作用） |
| `delta_max` | 0.6 | 最大耗散比 |
| `x0_km, y0_km, z0_km` | 600 / 0 / 350 | 泡中心 [km] |
| `sigma_perp_km` | 50.0 | 垂直 B 半宽 [km]（典型 30–100 km） |
| `sigma_para_km` | 550.0 | 沿 B 半长 [km]（典型 500–1000 km） |

### 5.7 OTH 雷达参数（`config.RADAR`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `mode` | `'monostatic'` | 工作模式 |
| `target_range_km` | 1169.0 | TX → 目标单程距离 [km] |
| `sigma_rcs_m2` | 5.0 | 目标 RCS [m²]（10 MHz 飞机谐振区） |
| `two_way` | True | 输出双程时延 |
| `sigma0_ground_dB` | -20.0 | 地面归一化 RCS [dBsm]（海：-20，陆：-25） |

### 5.8 射线追踪参数（`config.RT` / `config.RT_3D`）

| 键 | RT 默认值 | RT_3D 默认值 | 说明 |
|----|-----------|--------------|------|
| `ds_km` | 0.5 | 0.5 | RK4 标称步长 [km] |
| `ds_min_km` | 0.1 | 0.1 | 自适应步长最小值 [km] |
| `ds_max_km` | 2.0 | 2.0 | 自适应步长最大值 [km] |
| `beta_min/max` | 5 / 85 | 5 / 85 | 扇形扫射仰角范围 [°] |
| `n_fan` | 33 | 33 | 扇形扫射射线数 |
| `z_stop_km` | 550.0 | 550.0 | 射线上升截断高度 [km] |
| `max_steps` | 4000 | 4000 | 每条射线最大 RK4 步数 |
| `earth_flat` | — | True | Jones 1975 展平（仅 RT_3D） |

### 5.9 P2P 变分参数（`config.P2P` / `config.P2P_3D`）

| 键 | P2P 默认值 | P2P_3D 默认值 | 说明 |
|----|------------|---------------|------|
| `n_init` | 18 | **6** | 初始仰角数（P2P_3D 用 fan 热启动，取最近 RX 的前 6 条） |
| `n_ctrl` | 30 | — | 控制点数（P2P 仅） |
| `alpha_km` | 0.5 | — | 梯度下降步长 [km]（P2P 仅） |
| `k_spring` | 0.1 | — | 平滑弹簧系数（P2P 仅） |
| `max_iter` | 500 | **12** | 最大迭代次数（P2P_3D 发散时提前退出） |
| `tol_km` | 0.05 | 5.0 | 收敛准则 [km] |
| `d_beta` | — | 0.05 | Jacobian 仰角步 [°]（P2P_3D 仅） |
| `d_az` | — | 0.10 | Jacobian 方位步 [°]（P2P_3D 仅） |
| `jacobian` | — | `forward` | 三维 P2P Jacobian 差分方式：`forward` 或 `central` |
| `clust_h_km` | 10.0 | 10.0 | 去重高度容差 [km] |
| `clust_tau_ms` | 0.05 | 0.05 | 去重时延容差 [ms] |
| `n_workers` | 1 | — | 并行进程（0=自动；脚本需 `__main__`） |

### 5.10 流管参数（`config.TUBE_TRACER` / `config.TUBE_3D`）

| 键 | TUBE_TRACER | TUBE_3D | 说明 |
|----|-------------|---------|------|
| `delta_beta_deg` | 0.5 | 0.5 | 仰角步长 [°] |
| `delta_az_deg` | — | 0.5 | 方位扰动步长 [°]（三维流管） |
| `n_tube_rays` | 80 | 80 | 最大射线数 |
| `x_tgt_tol_km` | 80.0 | 80.0 | 目标落点搜索窗口 [km] |
| `T_pulse_ms` | 0.5 | 0.5 | 雷达脉冲宽度 [ms] |
| `L_cross_km` | 100.0 | 100.0 | 方位向波束足迹 [km] |
| `newton_tol_km` | 1.0 | — | Newton 精化收敛 [km] |
| `newton_max_iter` | 5 | — | Newton 最大迭代次数 |

### 5.11 散射参数（`config.SCATTER`，Phase 6）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `Cs_rel` | 0.01 | 相对 RMS 起伏幅度 dNe/Ne（无量纲） |
| `p` | 3.0 | 三维 von Karman 谱指数（F 层典型 2.5–4.0） |
| `L_outer_km` | 10.0 | 外尺度（相关长度）[km] |
| `L_inner_m` | 0.1 | 内尺度截断 [m] |
| `DeltaL_km` | 100.0 | 散射层路径厚度 [km] |
| `z_eff_km` | 300.0 | Rytov Fresnel 滤波器有效观测距离 [km] |

### 5.12 相位屏参数（`config.PHASE_SCREEN`，Phase 7）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `N_screen` | 10 | MPS 相位屏数量 |
| `n_pts` | 512 | 横向网格点数（n_pts × n_pts） |
| `dx_m` | 200.0 | 横向网格间距 [m] |
| `seed` | None | 随机种子（None=随机；设整数以复现） |

### 5.13 地磁参数（`config.GEOMAG`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `fH_MHz` | 1.197 | 回旋频率 [MHz]（IGRF-14，32.5N 120E 300km） |
| `dip_deg` | 48.7 | 地磁倾角 [°]（正下倾） |
| `decl_deg` | -5.5 | 地磁偏角 [°]（负=偏西） |
| `enable_OX` | False | O/X 磁离子分裂总开关 |

> 更新地磁参数：`GEOMAG.update(get_geomag(IRI_LAT, IRI_LON, alt_km=300, dt=IRI_DT))`

### 5.14 其他参数

**D 层吸收（`config.ABSORPTION`）**：

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | 开关 |
| `A0` | 500.0 | 吸收系数 [dB·MHz²]（Pederick & Cervera 2014） |
| `chi_deg` | 60.0 | 太阳天顶角 [°]（60=白天，90=夜间=0 dB） |

**Spread-F（`config.SPREAD_F`）**：

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | 开关 |
| `Cs` | 1e-3 | 谱强度 |
| `p` | 3.0 | 幂律谱指数 |
| `h_screen_km` | 300.0 | 相位屏中心高度 [km] |
| `L0_km` | 50.0 | 外尺度 [km] |

**PE/SSF（`config.PE`）**：

| 键 | 默认值 | 说明 |
|----|--------|------|
| `dx_km` | 0.5 | 传播步长 [km] |
| `dz_m` | λ/4 = 7.5 m | 垂直采样（@10 MHz） |
| `n_pml` | 60 | PML 吸收层厚度 [格点] |
| `sigma_pml` | 0.4 | PML 最大衰减系数 |
| `w0_km` | 20.0 | Gaussian 波束腰 [km] |
| `z_pe_min/max_km` | 200 / 500 | PE 域高度范围 [km] |
| `store_history` | False | 存储所有 x 截面（调试，~GB） |
| `earth_flat` | True | n_eff = n*(1+z/R_E) 修正 |

---

## 6. 输出格式

### 6.1 mode_results 字典键

| 键 | 单位 | 说明 | 可用路径 |
|----|------|------|---------|
| `label` | — | 模式标签（`1F_low`, `2F_O`, `Es`, …） | 全部 |
| `tau_ms` | ms | 单程群时延 | 全部 |
| `tau_2way_ms` | ms | 双程时延（雷达/流管） | 雷达/流管 |
| `delta_tau_ms` | ms | Es 或 PE 引入的时延扩展 | P2P |
| `tau_spread_ms` | ms | 管内时延扩展（相邻射线时延差） | 流管 |
| `Pr_W` | W | 接收功率（线性） | 全部 |
| `Pr_dBW` | dBW | 接收功率（对数） | 全部 |
| `h_reflect_km` | km | 反射点高度 | 全部 |
| `group_path_km` | km | 单程群路径 | 全部 |
| `geom_path_km` | km | 几何射线路径（区别于群路径） | 射线追踪 |
| `beta_deg` | ° | 出发仰角 | 全部 |
| `phi_deg` | ° | 链路方位角 | 全部 |
| `wave_mode` | — | `'O'`/`'X'`/`'iso'` | 全部 |
| `F_focus` | — | 聚焦因子 | 流管 |
| `F_focus_el` | — | 仰角平面聚焦因子 | 三维流管 |
| `F_focus_az` | — | 方位平面聚焦因子 | 三维流管 |
| `F_focus_3d` | — | 三维几何平均聚焦因子 | 三维流管 |
| `A_tube_km2` | km² | 流管地面投影面积 | 流管 |
| `abs_dB` | dB | D 层单程吸收诊断 | 吸收开启时/3D 功率链 |
| `abs_total_dB` | dB | D 层总吸收（通信 1x，雷达 2x） | 吸收开启时/3D 功率链 |
| `y_land_km` | km | 落地横向偏移 | 三维 |
| `azimuth_deflect_deg` | ° | 方位偏转角 arctan(y/x) | 三维 |
| `points` | km | 控制点数组 (n_ctrl+2,2) | P2P |

### 6.2 mode_summary 字典（`model.mode_summary`）

| 键 | 类型 | 说明 |
|----|------|------|
| `n_modes` | int | 传播模式总数 |
| `main_label` | str | 主模式标签（P-D 谱峰功率对应） |
| `main_tau_ms` | float | 主模式单程时延 [ms] |
| `main_Pr_dBW` | float | 主模式功率 [dBW] |
| `main_F_focus` | float | 主模式聚焦因子（P2P 路径时为 1.0） |
| `tau_spread_main` | float | 主模式时延扩展 [ms] |
| `mode_pairs` | list | O/X 配对列表 |

### 6.3 三维模式额外输出字段

三维 P2P 和三维流管返回的 ray dict 中额外包含：

```python
{
    'x_land_km'          : float,    # 落地 x [km]
    'y_land_km'          : float,    # 落地横向偏移 [km]（沿链路=x，横向=y）
    'azimuth_deflect_deg': float,    # arctan(y_land / x_land) [°]
    'beta_deg_tx'        : float,    # 发射仰角 [°]
    'az_deg_tx'          : float,    # 发射方位偏转 [°]（0=沿链路）
    'beta_deg_rx'        : float,    # 到达仰角 [°]
    'az_deg_rx'          : float,    # 到达方位偏转 [°]
    'h_reflect_km'       : float,    # 反射点物理高度（含展平修正）[km]
    'geom_path_km'        : float,    # 几何路径 integral(ds) [km]
    'group_path_km'       : float,    # 群路径 integral(ds/n) [km]
    'trajectory_3d'      : list,     # [(x,y,z,kx,ky,kz), ...] 全轨迹
}
```

### 6.4 Born/Rytov 额外输出（Phase 6）

三维主流程通过 `models/scatter_coupler.py` 将散射 regime 诊断挂到 mode dict 上，默认不直接改写 `Pr_W`：

```python
{
    'scatter_method'    : str,    # 'Born' / 'Rytov' / 'MPS_required'
    'sigma_phi_sq'      : float,  # 相位方差 [rad^2]
    'sigma_phi_rad'     : float,  # RMS 相位 [rad]
    'scatter_Ne0_m3'    : float,  # 沿三维轨迹估计的代表背景 Ne
    'scatter_Cs_rel'    : float,
    'scatter_L_outer_km': float,
    'scatter_DeltaL_km' : float,
    'mps_required'      : bool,
}
```

```python
# scatter_rytov.rytov_full() 返回
{
    'S4'            : float,   # 闪烁指数（0–1 弱散射，>0.6 强散射）
    'sigma_phi_rad' : float,   # RMS 相位抖动 [rad]
    'sigma_tau_ms'  : float,   # RMS 时延抖动 [ms]
    'fresnel_r_km'  : float,   # Fresnel 尺度 [km]
}

# scatter_born.born_sigma_v() 返回
sigma_v  # float, 体散射截面 [m^-1]
```

---

## 7. 场景运行示例

### 7.1 main.py 中的 9 个场景

| 函数 | 电离层 | 特殊功能 | 输出文件前缀 |
|------|--------|----------|-------------|
| `run_baseline()` | IRI 背景 | — | `ray_fan_baseline`, `pd_baseline` |
| `run_with_tid()` | IRI + TID | 多方向叠加 | `ray_fan_tid`, `pd_tid` |
| `run_with_es()` | IRI + Es | Hao 三段式 | `ray_fan_es`, `pd_es` |
| `run_with_bubble()` | IRI + 3-D 场向泡 | Gaussian 磁场对齐耗空 | `ray_fan_bubble`, `pd_bubble` |
| `run_full()` | IRI + TID + Es + 泡 | 全特效 | `ray_fan_full`, `pd_full` |
| `run_with_OX()` | IRI + TID | O/X 磁离子分裂 | `ray_fan_OX`, `pd_OX` |
| `run_radar_baseline()` | IRI | 雷达方程 | `ray_fan_radar`, `pd_radar` |
| `run_with_spreadf()` | IRI + TID + Spread-F | 幂律相位屏 | `ray_fan_spreadf`, `pd_spreadf` |
| `run_tube_radar()` | IRI | Coleman 流管 | `pd_tube_radar` |

### 7.2 代码示例

```python
import config as cfg
from models.hybrid_model import HybridPropagationModel

# ── 最简（IRI 背景，通信模式）────────────────────────────────────────────────
model = HybridPropagationModel(
    iono_params  = {'iri_params': {'dt': cfg.IRI_DT,
                                    'lat': cfg.IRI_LAT,
                                    'lon': cfg.IRI_LON}},
    radar_params = {'freq_MHz': cfg.FREQ_MHZ, 'Pt_W': cfg.PT_W,
                    'Gt': cfg.GT, 'Gr': cfg.GR},
)
modes, tau_ax, pd_W, main = model.compute(cfg.TX_POS, cfg.RX_POS)

# ── 三维 P2P 直接调用 ─────────────────────────────────────────────────────────
from models.ionosphere_model import IonosphereModel
from models.ray_tracer_3d import RefractiveIndex3D
from models.point_to_point_3d import find_all_rays_3d

iono = IonosphereModel(tid_params={**cfg.TID, 'enable': True})
Ne, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
n3d   = RefractiveIndex3D(Ne, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ,
                           tid_params={**cfg.TID, 'enable': True},
                           earth_flat=True)
modes3d = find_all_rays_3d(cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ,
                            p2p_params=cfg.P2P_3D)
for r in modes3d:
    print("{} tau={:.3f}ms h={:.0f}km az_defl={:.2f}deg".format(
        r['label'], r['tau_ms'], r['h_reflect_km'], r['azimuth_deflect_deg']))

# ── 三维场向泡：n2 对比 ───────────────────────────────────────────────────────
bub_on  = {**cfg.BUBBLE_3D, 'enable': True}
n3d_bub = RefractiveIndex3D(Ne, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ,
                              bubble_3d_params=bub_on, earth_flat=True)
x_c, z_phys = cfg.BUBBLE_3D['x0_km'], cfg.BUBBLE_3D['z0_km']
z_cart = z_phys + x_c**2 / (2 * cfg.RE_KM)
print("n2 off={:.4f} on={:.4f} delta={:.5f}".format(
    n3d.n2(x_c, 0, z_cart),
    n3d_bub.n2(x_c, 0, z_cart),
    n3d_bub.n2(x_c, 0, z_cart) - n3d.n2(x_c, 0, z_cart)))

# ── Born/Rytov 散射快速评估 ───────────────────────────────────────────────────
from models.scatter_selector import scatter_phase_variance, select_scatter_method
from models.scatter_born import born_sigma_v
from models.scatter_rytov import rytov_full

Ne0 = 5e10   # F 层电子密度 [m^-3]
sp  = scatter_phase_variance(Ne0, cfg.SCATTER['Cs_rel'],
                              cfg.SCATTER['DeltaL_km'],
                              cfg.SCATTER['L_outer_km'], cfg.FREQ_MHZ)
print("sigma_phi^2={:.3f}  method={}".format(sp, select_scatter_method(sp)))

ryt = rytov_full(cfg.FREQ_MHZ, cfg.SCATTER['Cs_rel'], cfg.SCATTER['p'],
                  cfg.SCATTER['L_outer_km'], Ne0,
                  cfg.SCATTER['z_eff_km'], cfg.SCATTER['DeltaL_km'])
print("S4={:.4f}  sigma_phi={:.3f}rad  sigma_tau={:.2e}ms".format(
    ryt['S4'], ryt['sigma_phi_rad'], ryt['sigma_tau_ms']))

# ── MPS 强散射（等离子体泡） ──────────────────────────────────────────────────
from models.scatter_phase_screen import MultiPhaseScreen
mps = MultiPhaseScreen(freq_MHz=cfg.FREQ_MHZ, Cs_rel=0.30,
                        p=cfg.SCATTER['p'], L_outer_km=cfg.SCATTER['L_outer_km'],
                        Ne_m3=Ne0, DeltaL_km=cfg.SCATTER['DeltaL_km'],
                        N_screen=cfg.PHASE_SCREEN['N_screen'],
                        n_pts=128, dx_m=200.0, seed=0)
res = mps.propagate()
print("MPS S4={:.4f}  Pr_norm={:.3f}".format(res['S4'], res['Pr_norm']))
```

---

## 8. 验证脚本汇总

全部应输出 `0 fail(s)`。快速验证顺序：

```bash
# 基础（Phase 1-4，按顺序）
conda run -n pytorch_cpu python tests/_verify_part0.py
conda run -n pytorch_cpu python tests/_verify_part1.py
conda run -n pytorch_cpu python tests/_verify_part3.py
conda run -n pytorch_cpu python tests/_verify_part6.py
conda run -n pytorch_cpu python tests/_verify_radar.py
conda run -n pytorch_cpu python tests/_verify_OX.py
conda run -n pytorch_cpu python tests/_verify_absorption.py
conda run -n pytorch_cpu python tests/_verify_absorption_3d.py
conda run -n pytorch_cpu python tests/_verify_group_path.py
conda run -n pytorch_cpu python tests/_verify_p2p_jacobian.py
conda run -n pytorch_cpu python tests/_verify_scatter_coupling.py
conda run -n pytorch_cpu python tests/_verify_spreadf.py
conda run -n pytorch_cpu python tests/_verify_tube.py

# 高级（Phase 5-7）
conda run -n pytorch_cpu python tests/_verify_phase5.py
conda run -n pytorch_cpu python tests/_verify_phase6.py
conda run -n pytorch_cpu python tests/_verify_phase7.py
```

**关键验证数值汇总**：

| 脚本 | 关键指标 | 期望值 |
|------|---------|--------|
| `_verify_part1.py` | IRI Ne 峰值，TID 振幅，泡耗散 | 定量断言通过 |
| `_verify_group_path.py` | `group_path=integral(ds/n)`，`tau=group_path/c` | 0 fail(s) |
| `_verify_absorption_3d.py` | 三维通信 1x / 雷达 2x D 层吸收 | 0 fail(s) |
| `_verify_scatter_coupling.py` | Born/Rytov/MPS regime 选择与字段输出 | 0 fail(s) |
| `_verify_p2p_jacobian.py` | forward vs central Jacobian 主模式一致性 | dtau < 0.2 ms |
| `_verify_tube.py` | 流管 tau vs P2P tau 差值 | < 0.05 ms |
| `_verify_tube.py` | F_focus 范围 | (0.1, 10.0) |
| `_verify_phase5.py` | 三维 TID az=30° 横向偏移 | > 0.05 km |
| `_verify_phase6.py` | Check 1: 三维射线 tau vs 二维差 | < 0.1 ms |
| `_verify_phase6.py` | Check 3: Born sigma_v 范围 | [1e-20, 1e-5] m^-1 |
| `_verify_phase6.py` | Check 5: Rytov S4 | [0, 0.6] |
| `_verify_phase7.py` | Check 1: 泡中心 delta_n2 | > 0.001 |
| `_verify_phase7.py` | Check 2: 三维流管模式数 | >= 1 |
| `_verify_phase7.py` | Check 3: BG sigma_v(asp=0) | [1e-30, 1.0] m^-1 |
| `_verify_phase7.py` | Check 4: Rino 屏 ratio | [0.1, 10.0] |
| `_verify_phase7.py` | Check 5: MPS S4 单调性 | S4(Cs=0.10) > S4(Cs=0.005) |

---

## 9. 常见修改方法

### 9.1 更换链路（TX/RX 位置和距离）

```python
# config.py 中修改
TX_LAT, TX_LON = 25.0, 105.0      # 新发射机位置
RX_RANGE = 800.0                   # 新目标距离
RX_POS   = (800.0, 0.0)
RADAR['target_range_km'] = 800.0
IRI_LAT  = 27.0                    # 新路径中点纬度
BG_X_MAX = 1000.0                  # 须大于 RX_RANGE + ~200
BG_X = np.arange(BG_X_MIN, BG_X_MAX + BG_DX, BG_DX)
# 更新 GEOMAG:
# from utils import get_geomag
# GEOMAG.update(get_geomag(IRI_LAT, IRI_LON, alt_km=300, dt=IRI_DT))
```

### 9.2 更换频率

```python
# config.py
FREQ_MHZ = 15.0
# PE['dz_m'] 随 _LAM_M/4 自动更新
# 若使用 O/X，更新 GEOMAG['fH_MHz']（参见 get_geomag）
```

### 9.3 TID 多方向叠加（Phase 5）

```python
# config.py 中设置两个方向的 MSTID：
TID = {
    **TID,
    'enable':           True,
    'n_components':     2,
    'az_deg_list':      [0.0, 45.0],       # 两个分量方向（deg 从正北）
    'amplitude_list':   [0.10, 0.05],
    'period_s_list':    [2400.0, 1800.0],
    'lambda_h_km_list': [300.0, 200.0],
}
```

### 9.4 开启三维场向等离子体泡

```python
# config.py
BUBBLE_3D = {**BUBBLE_3D, 'enable': True, 'delta_max': 0.8,
             'sigma_perp_km': 40.0, 'sigma_para_km': 600.0}

# 在 RefractiveIndex3D 中启用
n3d = RefractiveIndex3D(Ne, BG_X, BG_Z, FREQ_MHZ,
                         bubble_3d_params=BUBBLE_3D, earth_flat=True)
```

### 9.5 调整散射强度并分析方法选择

```python
from models.scatter_selector import scatter_phase_variance, select_scatter_method

for Cs in [0.001, 0.01, 0.05, 0.10, 0.30]:
    sp = scatter_phase_variance(5e10, Cs, 100.0, 10.0, 10.0)
    print("Cs={:.3f}: sigma_phi^2={:.2e}  -> {}".format(
        Cs, sp, select_scatter_method(sp)))
```

### 9.6 Booker-Gordon Es 方位角灵敏度分析

```python
from models.es_model import booker_gordon_sigma_v, es_scatter_valid
import numpy as np

L_para = 10000.0  # 10 km 场向尺度
for asp in [0.0, 0.01, 0.05, 0.1, 1.0, 5.0, 90.0]:
    sv = booker_gordon_sigma_v(10.0, 0.1, 2e10, 300.0, L_para, 1.5, asp)
    ok = es_scatter_valid(asp, L_para, 10.0)
    print("aspect={:6.2f}deg  sigma_v={:.2e}  valid={}".format(asp, sv, ok))
```

### 9.7 MPS 参数扫描（强散射分析）

```python
from models.scatter_phase_screen import MultiPhaseScreen

for Cs in [0.01, 0.05, 0.10, 0.20, 0.30]:
    mps = MultiPhaseScreen(freq_MHz=10.0, Cs_rel=Cs, p=3.0,
                            L_outer_km=10.0, Ne_m3=5e10,
                            DeltaL_km=100.0, N_screen=10,
                            n_pts=128, dx_m=200.0, seed=42)
    res = mps.propagate()
    print("Cs={:.2f}  S4={:.4f}  Pr_norm={:.3f}".format(
        Cs, res['S4'], res['Pr_norm']))
```

### 9.8 加速技巧

```python
# P2P 快速模式（减少初始路径数和迭代次数）
p2p_fast = {**cfg.P2P, 'n_init': 6, 'max_iter': 100}
modes, *_ = model.compute(p2p_params=p2p_fast)

# P2P 并行加速（Windows 必须有 __main__ 保护）
if __name__ == '__main__':
    modes, *_ = model.compute(p2p_params={**cfg.P2P, 'n_workers': 0})

# 减小 MPS 网格（测试用）
mps_fast = MultiPhaseScreen(..., n_pts=64, N_screen=5)

# 减少三维流管射线数
tube3d_fast = {**cfg.TUBE_3D, 'n_tube_rays': 40, 'delta_beta_deg': 1.0}
```

### 9.9 地磁参数更新

```python
# 针对新链路重算 IGRF-14 地磁参数（需 ppigrf 已安装）
from utils import get_geomag
import config as cfg

new_geomag = get_geomag(lat=cfg.IRI_LAT, lon=cfg.IRI_LON,
                         alt_km=300, dt=cfg.IRI_DT)
cfg.GEOMAG.update(new_geomag)
print(cfg.GEOMAG)
# {'fH_MHz': ..., 'dip_deg': ..., 'decl_deg': ..., 'enable_OX': False}
```

---

## 10. 性能参考

测试环境：Windows 11，单线程，pytorch_cpu conda 环境（_bench.py 实测）

### 10.1 实测基准（_bench.py）

```
射线扇面（33 条，max_steps=4000）    :  ~4.5 s
3-D P2P（n_init=6，fan 热启动）      :  ~9.4 s
场景 0 合计（扇面 + P2P）            : ~14.0 s
全部 9 个场景                        :  ~3–5 min
```

### 10.2 各操作耗时

| 操作 | 典型耗时 | 快速替代 |
|------|---------|---------|
| M1: IRI 背景 | ~1 s | — |
| M1: IRI + 烘焙 TID | ~1.5 s | — |
| M2b: 2-D P2P（n_init=18，iter=500） | 30–60 s | n_init=6, iter=100: ~5 s |
| M2d: **3-D P2P（n_init=6，fan 热启动）** | **~9 s** | — |
| M2d: 3-D P2P（n_init=18，无热启动，旧版） | 3–10 min | 勿用 |
| M2c: 2-D 流管（80 射线，Newton=True） | 10–30 s | n_tube_rays=40: ~5 s |
| A3b: 3-D 流管（80 射线 + az 扰动） | 20–60 s | n_tube_rays=40 |
| C1/C2: Born/Rytov（单点） | < 0.1 s | — |
| C5: MPS（n_pts=512，N_screen=10） | 5–20 s | n_pts=64, N_screen=5: ~0.5 s |
| M4: PE/SSF（100 km 域，dz=7.5 m） | 10–100 s | — |
| 完整 run_tube_radar() | ~30 s | — |
| 完整 run_full()（TID + Es + 泡） | ~25 s | — |

### 10.3 关键性能优化（已实现）

| 优化 | 加速比 | 说明 |
|------|--------|------|
| `_bilin2` 快速双线性插值 | ~15x | 取代 scipy RGI，O(1) 直接数组索引，避免 scipy 开销 |
| 预计算 Cartesian n² 梯度 | ~7x | `__init__` 时计算 `_dn2_dx_arr/_dn2_dz_arr`，`grad_n2()` 直接查表 |
| 烘焙 TID（baked TID） | ~5x | TID 预合并入 Ne_2d，`RefractiveIndex3D` 不传 `tid_params`，避免每步 4 次 RGI |
| 单侧 Jacobian | ~1.4x | 3 条射线/迭代（vs 中心差分 5 条）：P2P_3D 节省 40% 射线追踪 |
| 早退出 | ~2x | 第 3 次迭代后误差不减则放弃该种子 |
| Fan 热启动 | 总计 **~23x** | P2P 从 221 s 降至 9.4 s；fan 按 `|x_land - x_rx|` 排序取最近 n_init 条为种子 |

> **烘焙 TID 用法**（main.py 中 TID 场景的标准做法）：
> ```python
> # 正确：TID 预合并进 Ne_2d，3-D RT 不感知 TID，走快速路径
> iono  = IonosphereModel(tid_params={**cfg.TID_3D, 'enable': True})
> Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
> n3d   = RefractiveIndex3D(Ne_2d, ...)          # 不传 tid_params
>
> # 错误（卡顿）：每步 n2() 调用都做 4 次 RGI，33 条射线 × 4000 步 = 5.3 万次 RGI
> n3d   = RefractiveIndex3D(Ne_2d, ..., tid_params={**cfg.TID_3D, 'enable': True})
> ```

---

## 11. 参考文献

| # | 文献 | 对应模块 |
|---|------|---------|
| [1] | Hao, Y. et al. (2017). Sporadic-E layer scintillation. *Radio Science* | M3 Es 三段式 |
| [2] | Hooke, W.H. (1968). Ionospheric irregularities produced by internal atmospheric gravity waves. *J. Atm. Terr. Phys.*, 30, 795–823 | M1 TID |
| [3] | Carrano, C. et al. (2020). Wide-angle PE propagation. *Radio Science* | M4 PE/SSF |
| [4] | Nosikov, I.A. et al. (2020). Variational P2P solver. *Radio Science*, 55 | M2b P2P |
| [5] | Budden, K.G. (1961). *Radio Waves in the Ionosphere*. Cambridge | M2a Appleton-Hartree |
| [6] | Pederick, L.H. & Cervera, M.A. (2014). HF absorption. *Radio Science*, 49, 81–93 | D 层吸收 |
| [7] | Rino, C.L. (1979). Power law phase screen. *Radio Science*, 14(6), 1135–1145 (Parts 1 & 2) | C3 Rino 屏 |
| [8] | Jones, R.M. & Stephenson, J.J. (1975). A versatile 3-D ray tracing program. OT Report 75-76 | M2d 地球展平 |
| [9] | Cervera, M.A. & Harris, T.J. (2014). Modeling ionospheric disturbance features. *JGR*, 119, 431–440 | M2d Haselgrove-3D |
| [10] | Coleman, C.J. (1997). A ray tracing formulation for OTH radar. *Radio Science*, 32(1) | A3 流管 |
| [11] | Coleman, C.J. (1998). A model of HF sky wave radar clutter. *Radio Science*, 33(4) | A3 流管功率 |
| [12] | Galushko, V.G. et al. (2013). Born approximation in HF scattering. *Radio Science*, 48, 577–586 | C1 Born |
| [13] | Yeh, K.C. & Liu, C.H. (1982). Radio wave scintillations in the ionosphere. *Proc. IEEE*, 70(4) | C2 Rytov |
| [14] | Knepp, D.L. (1983). Multiple phase-screen calculation. *Proc. IEEE*, 71(6), 722–737 | C5 MPS |
| [15] | Booker, H.G. (1956). Theory of scattering by nonisotropic irregularities (aurora). *J. Atm. Terr. Phys.*, 8, 204–221 | C4 Booker-Gordon |
| [16] | Booker, H.G. & Gordon, W.E. (1950). Theory of radio scattering in troposphere. *Proc. IRE*, 38 | C4 BG 源头 |
| [17] | Ma, G. et al. (2023). Three-dimensional ray tracing through equatorial plasma bubbles. *GPS Solutions*, 27, 45 | B4 场向泡 |
| [18] | Martinis, C. (2023). Equatorial plasma bubbles: 3-D modeling. *JGR*, 128, e2023JA031820 | B4 场向泡 |
| [19] | Tatarski, V.I. (1961). *Wave Propagation in a Turbulent Medium*. McGraw-Hill | C1/C2 湍流理论 |
| [20] | Jiang, C. et al. (2020). Plasma bubble modeling. | M1 等离子体泡 |
