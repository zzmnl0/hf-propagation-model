# HF Shortwave Hybrid Propagation Model

研究项目 李慧敏-22：复杂电离层环境短波反射/散射混合传播机理及建模研究

**测试链路**：TX @ 30°N 120°E，目标距离 1169 km，f = 10 MHz  
**代码仓库**：`https://github.com/zzmnl0/hf-propagation-model`  
**最新提交**：Phase 3 全部完成（D 层吸收 + 扩展 F 相位屏 + 实测比对框架）

---

## 目录

1. [快速上手](#1-快速上手)
2. [目录结构](#2-目录结构)
3. [整体数据流](#3-整体数据流)
4. [模块详解](#4-模块详解)
5. [参数速查表](#5-参数速查表)
6. [输出格式](#6-输出格式)
7. [场景运行示例](#7-场景运行示例)
8. [验证脚本](#8-验证脚本)
9. [常见修改方法](#9-常见修改方法)
10. [性能参考](#10-性能参考)
11. [参考文献](#11-参考文献)

---

## 1. 快速上手

```bash
# 完整运行所有场景（输出到 output/）
conda run -n pytorch_cpu python main.py

# 单独运行某一场景（在 main.py 底部修改调用列表，或直接 import）
conda run -n pytorch_cpu python -c "
import sys, os; sys.path.insert(0, '.')
from main import run_baseline
run_baseline()
"

# 验证所有 Phase
conda run -n pytorch_cpu python tests/_verify_part1.py
conda run -n pytorch_cpu python tests/_verify_OX.py
conda run -n pytorch_cpu python tests/_verify_radar.py
conda run -n pytorch_cpu python tests/_verify_absorption.py
conda run -n pytorch_cpu python tests/_verify_spreadf.py

# 背景电子密度图
conda run -n pytorch_cpu python viz/plot_ne_background.py
```

**依赖包**：`numpy` `scipy` `matplotlib` `iri2016` `ppigrf`

> Windows 终端为 GBK 编码，`print()` 中只能使用 ASCII 字符，中文注释不影响运行。

---

## 2. 目录结构

```
code/
├── config.py                    # 所有物理参数（唯一全局配置入口）
├── utils.py                     # 物理工具函数（EM、功率、坐标、D层吸收）
├── main.py                      # 场景入口（7个场景函数 + __main__）
│
├── models/
│   ├── ionosphere_model.py      # M1: IRI + TID / Es / 等离子体泡 / Spread-F
│   ├── ray_tracer.py            # M2a: Haselgrove RK4 + RefractiveIndexAH (O/X)
│   ├── point_to_point.py        # M2b: Nosikov 2020 变分 P2P 求解器
│   ├── es_model.py              # M3: Hao 2017 Es 三段式模型
│   ├── pe_propagator.py         # M4: Carrano 2020 PE/SSF 等离子体泡散射
│   ├── spread_f_model.py        # M3b: Rino 1979 幂律相位屏（扩展 F）
│   └── hybrid_model.py          # M5-M7: 完整 pipeline（P-D 谱 + 主模式识别）
│
├── analysis/
│   └── compare_pd.py            # 实测 vs 模型 P-D 谱比对框架（CLI 工具）
│
├── viz/
│   ├── plot_utils.py            # 共享绘图函数（射线扇、P-D 谱、CSV 导出）
│   └── plot_ne_background.py    # 4 面板 Ne 背景图
│
├── tests/
│   ├── _verify_part0.py         # 导入 / 配置健全性
│   ├── _verify_part1.py         # IonosphereModel 完整验证
│   ├── _verify_part3.py         # P2P 变分求解验证
│   ├── _verify_part4.py         # Es 模型验证
│   ├── _verify_part5.py         # PE/SSF 验证
│   ├── _verify_part6.py         # 混合模型 pipeline 验证
│   ├── _verify_radar.py         # Phase 1: OTH 雷达双程几何验证
│   ├── _verify_OX.py            # Phase 2: O/X 磁离子分裂验证
│   ├── _verify_absorption.py    # Phase 3: D 层吸收验证
│   └── _verify_spreadf.py       # Phase 3: 扩展 F 相位屏验证
│
└── output/                      # 生成 PNG / CSV（gitignored）
```

---

## 3. 整体数据流

```
config.py  （所有参数集中于此）
    |
    v
IonosphereModel.build_Ne_field(x, z, t)
    |-- 1. IRI-2016 一维背景 profile     [iri2016]
    |-- 2. TID 扰动  (Hooke 1968)        [enable=True]
    |-- 3. Es 薄层   (Hao 2017)          [enable=True]
    |-- 4. 等离子体泡 (Gaussian 耗散)    [enable=True]
    |-- 5. Spread-F  (Rino 1979 相位屏)  [enable=True]
    |
    v  Ne_2d (Nx, Nz) [m^-3]
    |
    v
RefractiveIndex  (各向同性) 或  RefractiveIndexAH  (Appleton-Hartree O/X)
    |
    v
find_all_rays_p2p(tx, rx, n_model, freq,
                  wave_mode='both'|None, geomag=GEOMAG)
    |-- wave_mode=None  -> 各向同性，wave_mode 字段='iso'
    |-- wave_mode='both' -> 分别用 O / X 折射率跑变分，
    |                        label 追加 _O / _X 后缀
    |-- 对每个模式：梯度下降变分 + 去重 + classify_mode
    |
    v  rays: list[dict]  (points, tau_ms, h_reflect_km, wave_mode, label, ...)
    |
    v
对每条 ray：
    |-- 功率初值:
    |     comm 模式: Friis  Pr = Pt*Gt*Gr / L_free
    |     radar 模式: 雷达方程  Pr = Pt*Gt*Gr*lambda^2*sigma / ((4pi)^3 * R^4)
    |
    |-- D 层吸收修正 [absorption.enable=True]:
    |     A_dB = A0*cos(chi)^0.75 / ((f+fH_L)^2 * sin(beta))
    |     radar: 2 次穿越 D 层（双程）
    |
    |-- Es 修正  [EsLayerModel, es.enable=True]:
    |     -> Pr_W 修正, delta_tau_ms 累加
    |
    |-- 等离子体泡 PE/SSF  [bubble.enable=True]:
    |     -> delta_tau_ms 累加, Pr_W *= P_out/P_in
    |
    v  mode_results: list[dict]  (含 tau_2way_ms, wave_mode, phi_deg 等新字段)
    |
    v
build_pd_spectrum -> identify_main_mode
    -> (mode_results, tau_axis, pd_W, main_mode)
```

---

## 4. 模块详解

### M1 — IonosphereModel

**文件**：`models/ionosphere_model.py`  
**类**：`IonosphereModel`

```python
iono = IonosphereModel(
    iri_params      = {'dt': IRI_DT, 'lat': IRI_LAT, 'lon': IRI_LON},
    tid_params      = {**TID,      'enable': True},   # 可选
    es_params       = {**ES,       'enable': True},   # 可选
    bubble_params   = {**BUBBLE,   'enable': True},   # 可选
    spread_f_params = {**SPREAD_F, 'enable': True},   # 可选
    freq_MHz        = 10.0,
)
Ne_2d, n_2d = iono.build_Ne_field(x_array, z_array, t=0.0)
```

**叠加顺序**（各层均可独立开关）：

| 步骤 | 扰动 | 物理模型 | 参数组 |
|------|------|----------|--------|
| 1 | IRI 背景 | iri2016.IRI() 一维 profile | `iri_params` |
| 2 | TID | Hooke (1968) 行进式电离层扰动 | `TID` |
| 3 | Es 薄层 | Hao (2017) 密度剖面 | `ES` |
| 4 | 等离子体泡 | Gaussian 耗散 (delta_max) | `BUBBLE` |
| 5 | Spread-F | Rino (1979) 幂律相位屏 | `SPREAD_F` |

---

### M2a — RefractiveIndex / RefractiveIndexAH

**文件**：`models/ray_tracer.py`

**类 `RefractiveIndex`**（各向同性）：

| 方法 | 说明 |
|------|------|
| `n(x, z)` | 单点折射率（标量） |
| `grad_n2(x, z)` | 中心差分梯度 (dn²/dx, dn²/dz) |
| `n_batch(pts)` | 矢量化，`pts` 形状 (M,2) -> (M,) |

**类 `RefractiveIndexAH`**（Appleton-Hartree O/X 波，继承 RefractiveIndex）：

```python
n_ah = RefractiveIndexAH(Ne_2d, x_km, z_km, freq_MHz,
                          wave_mode='O',   # 'O' | 'X'
                          geomag={'fH_MHz': 1.197, 'dip_deg': 48.7})
```

**Appleton-Hartree 公式**（Budden 1961 完整形式）：

```
n^2_{O,X} = 1 - X / (1 - YT^2/[2(1-X)] ∓ sqrt(YT^4/[4(1-X)^2] + YL^2))

X   = fp^2/f^2         (等离子频率参数)
Y   = fH/f             (回旋频率参数)
YT  = Y * sin(alpha)   (横向分量)
YL  = Y * cos(alpha)   (纵向分量)
alpha = pi/2 - dip     (射线与 B 的夹角，固定磁倾角近似)
O 波取 - (上符号)，X 波取 + (下符号)
```

**函数**：

| 函数 | 说明 |
|------|------|
| `trace_single_ray(tx, beta_deg, n_model, freq_MHz)` | Haselgrove RK4，返回 ray dict |
| `shoot_rays_fan(tx_pos, n_model, rt_params)` | 扇形扫角，返回 ray list |

---

### M2b — P2P 变分求解

**文件**：`models/point_to_point.py`

**主函数**：

```python
rays = find_all_rays_p2p(
    tx_km, rx_km, n_model, freq_MHz,
    p2p_params = P2P,
    wave_mode  = 'both',   # None=各向同性 | 'O' | 'X' | 'both'
    geomag     = GEOMAG,   # wave_mode != None 时必须提供
)
```

**算法**（Nosikov 2020）：

- 光程泛函 `S = sum n(mid_i) * |seg_i|`（离散 Fermat 积分）
- 高角射线：梯度下降 `pts -= alpha * grad_perp`
- 低角射线：符号翻转（鞍点 -> 稳定点）
- 弹簧力 `k_spring * (pts[i+1]-2*pts[i]+pts[i-1])` 保持平滑
- 去重：h 差 < `clust_h_km` 且 tau 差 < `clust_tau_ms`

**模式标签规则**：

| h_reflect | 条件 | 基础标签 | O/X 开启后 |
|-----------|------|----------|-----------|
| < 140 km | — | `Es` | `Es_O` / `Es_X` |
| 140–200 km | — | `E` | `E_O` / `E_X` |
| 200–300 km | tau < 5 ms | `1F_low` | `1F_low_O` / `1F_low_X` |
| 200–300 km | tau >= 5 ms | `1F_high` | `1F_high_O` / ... |
| >= 300 km | — | `2F` | `2F_O` / `2F_X` |

各向同性时 wave_mode 字段 = `'iso'`，标签无后缀。

**辅助函数**：

| 函数 | 说明 |
|------|------|
| `classify_mode(ray_dict)` | 返回含 O/X 后缀的标签 |
| `extract_es_params(points, h_Es_km)` | 射线与 Es 层的交叉参数 |
| `extract_bubble_entry(points, z_bot_km)` | 射线进入等离子体泡的参数 |

---

### M3 — Es 层模型

**文件**：`models/es_model.py`  **类**：`EsLayerModel`

**三段式分区**（foEs/f 为基准，Hao 2017）：

| foEs/f | 模式 | 处理 |
|--------|------|------|
| > 0.25 (`fr`) | reflect | 反射系数 ρ²，label='Es_reflect' |
| 0.10–0.25 | mixed | 线性插值，label='Es_mixed' |
| < 0.10 (`fs`) | scatter | 散射截面 σ₅，label='Es_scatter' |

---

### M3b — SpreadFModel（扩展 F）

**文件**：`models/spread_f_model.py`  
**参考**：Rino (1979) + Ding et al. (2021)

```python
from models.spread_f_model import SpreadFModel

sfm = SpreadFModel(Cs=1e-3, p=3.0, h_screen_km=300.0, L0_km=50.0, seed=42)
Ne_2d_new = sfm.apply(Ne_2d, x_km, z_km)
```

**物理模型**：

```
dNe(x, z) = Cs * phi(x) * Ne_bg(z) * G(z - h_screen)

phi(x) : 单位 RMS 幂律随机场，功率谱 S(k) ~ (k^2+k0^2)^(-(p+1)/2)
G(z)   : Gaussian 垂直包络（半宽 20 km）
k0     : 外尺度波数 = 2*pi / L0_km
```

验证结果：
- `Cs=0` → dNe 精确为零
- FFT 功率谱斜率 ≈ -(p+1) = -4（实测 -3.21，差值 < 1.0）
- dNe 随 Cs 单调增大

---

### M4 — PE/SSF 传播器

**文件**：`models/pe_propagator.py`  **类**：`PEPropagator`

**SSF 单步**（Carrano 2020 宽角 PE）：
- Step A（折射，空域）：`u_half = u * exp(i k0 (n-1) dx)`
- Step B（衍射，谱域）：FFT -> `kx_eff = sqrt(k0^2-kz^2) - k0` -> IFFT

**地球展平修正**（`PE['earth_flat']=True`）：
```
n_eff(x, z) = n(x, z) * (1 + z/R_E)
```
在 350 km 高度约修正 5.5%，关闭方式：`config.PE['earth_flat'] = False`

---

### M5-M7 — HybridPropagationModel

**文件**：`models/hybrid_model.py`  **类**：`HybridPropagationModel`

```python
model = HybridPropagationModel(
    iono_params       = {...},      # 见下方 iono_params 结构
    radar_params      = {...},      # freq_MHz, Pt_W, Gt, Gr
    radar_mode        = False,      # True = 雷达方程 + tau_2way_ms
    geomag_params     = None,       # None=各向同性; {**GEOMAG,'enable_OX':True}=O/X
    absorption_params = None,       # None=不修正; {enable:True,A0:500,...}=D层吸收
)
modes, tau_ax, pd_W, main = model.compute(tx_km, rx_km, t=0.0, p2p_params=P2P)
```

**iono_params 完整结构**：

```python
iono_params = {
    'iri_params':      {'dt': IRI_DT, 'lat': IRI_LAT, 'lon': IRI_LON},
    'tid_params':      {**TID,      'enable': True/False},
    'es_params':       {**ES,       'enable': True/False},
    'bubble_params':   {**BUBBLE,   'enable': True/False},
    'spread_f_params': {**SPREAD_F, 'enable': True/False},
}
```

**功率计算逻辑**：

```python
# 通信模式（radar_mode=False）
Pr_W = Pt * Gt * Gr * 10^(-L_free_dB/10)

# OTH 雷达模式（radar_mode=True）
Pr_W = Pt*Gt*Gr * lambda^2*sigma / ((4*pi)^3 * R^4)   # R = group_path_km（单程）
tau_2way_ms = 2 * tau_ms                                 # 双程时延
```

**D 层吸收（absorption_params.enable=True）**：

```python
# D 层吸收（Pederick & Cervera 2014）
A_dB = A0 * cos(chi)^0.75 / ((f + fH_L)^2 * sin(beta))
fH_L = GEOMAG['fH_MHz'] * cos(GEOMAG['dip_deg'])   # 纵向回旋频率

# 通信模式：1 次 D 层穿越
Pr_W *= 10^(-A_dB / 10)

# 雷达模式：双程 2 次穿越
Pr_W *= 10^(-2*A_dB / 10)
```

**模块级函数**：

| 函数 | 说明 |
|------|------|
| `build_pd_spectrum(mode_results)` | 高斯展宽叠加 P-D 谱 |
| `identify_main_mode(tau_axis, pd_W, mode_results)` | 峰值检测 + 映射最近模式 |

---

### 辅助工具

**`utils.py`** 中的主要函数：

| 函数 | 说明 |
|------|------|
| `free_space_loss_dB(D_km, freq_MHz)` | Friis 自由空间损耗 [dB] |
| `radar_equation_W(Pt, Gt, Gr, freq, gp_km, sigma)` | 单基雷达接收功率 [W] |
| `d_layer_absorption_dB(freq, beta, chi=60, A0=500, fH_L=0.789)` | D 层非偏吸收 [dB] |
| `get_geomag(lat, lon, alt_km, dt)` | ppigrf IGRF-14 地磁参数获取 |
| `haversine_km(lat1,lon1,lat2,lon2)` | 大圆距离 [km] |
| `bearing_deg(lat1,lon1,lat2,lon2)` | 初始方位角 [deg] |

**`analysis/compare_pd.py`** 实测比对工具（CLI）：

```bash
python analysis/compare_pd.py \
    --measured  data/pd_measured_20240115.csv \
    --scenario  baseline \
    --tau_tol   0.5 \
    --output    output/compare_baseline.png
```

测量数据 CSV 格式（必需列：`tau_ms`, `power_dBW`）：
```
tau_ms, power_dBW, freq_MHz, datetime, bearing_deg
3.90, -85.2, 10.0, 2024-01-15T04:00:00Z, 0.0
```

函数接口：

| 函数 | 说明 |
|------|------|
| `load_measured_pd(csv_path)` | 返回 (tau_ms, pd_dBW, meta) |
| `overlay_plot(...)` | 双面板比对图（P-D 谱叠加 + 模式标记） |
| `score_mode_match(modes_model, peaks_meas, tau_tol=0.5)` | 命中率、tau 误差、Pr 误差 |

---

## 5. 参数速查表

所有参数集中于 `config.py`，修改该文件即可全局生效。

### 测试链路

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `TX_LAT` | 30.0 | 发射机纬度 [°N] |
| `TX_LON` | 120.0 | 发射机经度 [°E] |
| `RX_RANGE` | 1169.0 | TX-目标水平距离 [km] |
| `FREQ_MHZ` | 10.0 | 工作频率 [MHz] |
| `PT_W` | 1000.0 | 发射功率 [W] |
| `GT`, `GR` | 1.0 | TX/RX 增益（线性，各向同性） |
| `IRI_DT` | 2020-06-01 12:00 | IRI 日期时间（正午，中等太阳活动） |
| `IRI_LAT` | 32.5°N | IRI 路径中点纬度 |
| `IRI_LON` | 120.0°E | IRI 路径中点经度 |

### OTH 雷达参数（`config.RADAR`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `mode` | `'monostatic'` | 单/双基工作模式 |
| `target_range_km` | 1169.0 | TX → 目标单程距离 [km] |
| `sigma_rcs_m2` | 5.0 | 目标 RCS [m²]（10 MHz 飞机谐振区典型值） |
| `two_way` | True | 输出双程时延 tau_2way_ms |
| `LINK_BEARING_DEG` | 0.0 | TX → 目标方位角 [°]（正北=0，正东=90） |

### 地磁参数（`config.GEOMAG`）

由 ppigrf（IGRF-14）在路径中点自动计算（32.5°N, 120°E, 300 km, 2020-01-01）：

| 键 | 默认值 | 说明 |
|----|--------|------|
| `fH_MHz` | 1.197 | 回旋频率 [MHz]（B=42766 nT） |
| `dip_deg` | 48.7 | 地磁倾角 [°] |
| `decl_deg` | -5.5 | 地磁偏角 [°]（负=偏西） |
| `enable_OX` | False | O/X 磁离子分裂总开关 |

更新方式：
```python
from utils import get_geomag
GEOMAG.update(get_geomag(IRI_LAT, IRI_LON, dt=IRI_DT))
```

### D 层吸收参数（`config.ABSORPTION`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | D 层吸收修正开关 |
| `A0` | 500.0 | 吸收系数 [dB·MHz²]（Pederick & Cervera 2014） |
| `chi_deg` | 60.0 | 太阳天顶角 [°]（60=白天典型中纬度，90=夜间=0 dB） |

> fH_L 不在此字典中，由 `hybrid_model.py` 从 GEOMAG 动态计算。

### 扩展 F 参数（`config.SPREAD_F`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | 扩展 F 相位屏开关 |
| `Cs` | 1e-3 | 相位谱强度（Rino 1979） |
| `p` | 3.0 | 幂律谱指数（典型 2.5–4.0） |
| `h_screen_km` | 300.0 | 相位屏中心高度 [km] |
| `L0_km` | 50.0 | 外尺度 [km] |

### TID 参数（`config.TID`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | TID 开关 |
| `lambda_h_km` | 300 | 水平波长 [km] |
| `T_s` | 2400 | 周期 [s]（40 分钟 MSTID） |
| `amplitude` | 0.10 | 峰值 dNe/Ne₀（0–1） |
| `I_dip_deg` | 50.0 | 地磁倾角 [°] |
| `H_km` | 60.0 | Chapman 标高 [km] |

### Es 层参数（`config.ES`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | Es 层开关 |
| `foEs_MHz` | 5.0 | Es 等离子频率 [MHz] |
| `h_Es_km` | 110.0 | Es 中心高度 [km] |
| `delta_h_m` | 115.0 | 半厚度 [m]（Hao 2017 最佳拟合） |
| `fr` | 0.25 | 反射阈值 foEs/f |
| `fs` | 0.10 | 散射阈值 foEs/f |

### 等离子体泡参数（`config.BUBBLE`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | 等离子体泡开关 |
| `delta_max` | 0.6 | 最大耗散比（0–1） |
| `x0_km` | 600.0 | 泡中心水平位置 [km] |
| `z0_km` | 350.0 | 泡中心高度 [km] |
| `Lx_km` | 100.0 | 水平半宽 [km] |
| `Lz_km` | 150.0 | 垂直半高 [km] |

### P2P 变分参数（`config.P2P`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `n_init` | 18 | 初始仰角扫描数 |
| `n_ctrl` | 30 | 路径控制点数 |
| `alpha_km` | 0.5 | 梯度下降步长 [km] |
| `k_spring` | 0.1 | 平滑弹簧系数 |
| `max_iter` | 500 | 最大迭代次数 |
| `tol_km` | 0.05 | 收敛准则 [km] |
| `clust_h_km` | 10.0 | 去重高度容差 [km] |
| `clust_tau_ms` | 0.05 | 去重时延容差 [ms] |
| `n_workers` | 1 | 并行进程数（0=自动，1=串行） |

### PE/SSF 参数（`config.PE`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `dx_km` | 0.5 | 传播步长 [km] |
| `dz_m` | lambda/4=7.5 m | 垂直采样（@10 MHz） |
| `n_pml` | 60 | PML 层厚度 [格点] |
| `earth_flat` | True | 地球展平修正 n_eff=n*(1+z/R_E) |
| `store_history` | False | 存储所有 x 截面（调试用） |

---

## 6. 输出格式

### mode_results 字典键

| 键 | 单位 | 说明 |
|----|------|------|
| `label` | — | 模式标签（见分类规则） |
| `tau_ms` | ms | 单程群时延 |
| `tau_2way_ms` | ms | 双程时延（雷达模式有值，通信模式为 None） |
| `delta_tau_ms` | ms | 时延扩展（来自 Es/PE 散射模型） |
| `Pr_W` | W | 接收功率（线性） |
| `Pr_dBW` | dBW | 接收功率（对数） |
| `h_reflect_km` | km | 反射高度（射线最高点） |
| `group_path_km` | km | 单程群路径长度 |
| `beta_deg` | deg | 出发仰角 |
| `phi_deg` | deg | 链路方位角（来自 LINK_BEARING_DEG） |
| `wave_mode` | — | `'O'` / `'X'` / `'iso'` |
| `points` | km | 变分控制点数组 (n_ctrl+2, 2)，用于路径可视化 |

### compute() 返回值

```python
modes, tau_ax, pd_W, main = model.compute(tx_km, rx_km)
# modes   : list[dict]         各传播模式
# tau_ax  : np.ndarray [ms]    P-D 谱横轴
# pd_W    : np.ndarray [W]     P-D 谱功率（高斯叠加）
# main    : dict | None        P-D 谱主峰对应的模式
```

### CSV 输出列

`save_modes_csv(scenario, modes, tau_ax, pd_W, freq_MHz, out_dir)` 输出：

```
scenario, label, freq_MHz, tau_ms, tau_2way_ms, delta_tau_ms,
Pr_W, Pr_dBW, pd_power_W, h_reflect_km, group_path_km, beta_deg, phi_deg
```

---

## 7. 场景运行示例

### 当前 7 个场景（main.py）

| 函数 | 电离层 | 特殊功能 | 输出文件前缀 |
|------|--------|----------|-------------|
| `run_baseline()` | IRI 仅背景 | — | `ray_fan_baseline`, `pd_baseline` |
| `run_with_tid()` | IRI + TID | — | `ray_fan_tid`, `pd_tid` |
| `run_with_es()` | IRI + Es | Es 三段式修正 | `ray_fan_es`, `pd_es` |
| `run_with_bubble()` | IRI + 等离子体泡 | PE/SSF 散射 | `ray_fan_bubble`, `pd_bubble` |
| `run_full()` | IRI + TID + Es + 泡 | 全特效 | `ray_fan_full`, `pd_full` |
| `run_with_OX()` | IRI + TID | O/X 磁离子分裂 | `ray_fan_OX`, `pd_OX` |
| `run_radar_baseline()` | IRI | OTH 雷达方程 | `ray_fan_radar`, `pd_radar` |
| `run_with_spreadf()` | IRI + TID + Spread-F | 幂律相位屏 | `ray_fan_spreadf`, `pd_spreadf` |

### 快速代码示例

```python
import config as cfg
from models.hybrid_model import HybridPropagationModel

# ── 最简运行（IRI 仅背景，通信模式）──────────────────────────────────────────
model = HybridPropagationModel(
    iono_params  = {'iri_params': {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON}},
    radar_params = {'freq_MHz': 10.0, 'Pt_W': 1000.0, 'Gt': 1.0, 'Gr': 1.0},
)
modes, tau_ax, pd_W, main = model.compute()

# ── OTH 雷达模式 ─────────────────────────────────────────────────────────────
model = HybridPropagationModel(
    iono_params  = {'iri_params': {...}},
    radar_params = {'freq_MHz': 10.0, 'Pt_W': 1e6, 'Gt': 1.0, 'Gr': 1.0,
                    'sigma_rcs_m2': cfg.RADAR['sigma_rcs_m2']},
    radar_mode   = True,
)
# mode['tau_2way_ms'] = 2 * mode['tau_ms']

# ── O/X 磁离子分裂 ────────────────────────────────────────────────────────────
model = HybridPropagationModel(
    iono_params   = {'iri_params': {...}, 'tid_params': {**cfg.TID, 'enable': True}},
    radar_params  = {'freq_MHz': 10.0, 'Pt_W': 1000.0, 'Gt': 1.0, 'Gr': 1.0},
    geomag_params = {**cfg.GEOMAG, 'enable_OX': True},
)
# 每个各向同性模式分裂为 *_O 和 *_X 两条，wave_mode 字段区分

# ── 开启 D 层吸收修正 ─────────────────────────────────────────────────────────
model = HybridPropagationModel(
    iono_params       = {...},
    radar_params      = {...},
    absorption_params = {**cfg.ABSORPTION, 'enable': True},
)
# 雷达模式自动应用 2x 吸收（双程两次穿越 D 层）

# ── 开启 Spread-F ────────────────────────────────────────────────────────────
iono_params = {
    'iri_params':      {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
    'tid_params':      {**cfg.TID,      'enable': True},
    'spread_f_params': {**cfg.SPREAD_F, 'enable': True, 'Cs': 5e-3},
}
model = HybridPropagationModel(iono_params, radar_params)

# ── 快速测试（减少 P2P 迭代）─────────────────────────────────────────────────
p2p_fast = {**cfg.P2P, 'n_init': 6, 'max_iter': 100}
modes, *_ = model.compute(p2p_params=p2p_fast)

# ── 开启 P2P 并行加速 ─────────────────────────────────────────────────────────
if __name__ == '__main__':
    p2p_par = {**cfg.P2P, 'n_workers': 0}   # 0=自动检测 CPU 数
    modes, *_ = model.compute(p2p_params=p2p_par)
```

---

## 8. 验证脚本

| 脚本 | 阶段 | 验证内容 | 关键数值 |
|------|------|----------|---------|
| `_verify_part1.py` | 基础 | IonosphereModel：Ne 值域、TID 振幅、Es 峰、泡耗散 | 定量断言 |
| `_verify_part3.py` | 基础 | 光程泛函、梯度、变分收敛、P2P 模式数 | 定量断言 |
| `_verify_part4.py` | 基础 | Es 分类（reflect/mixed/scatter）、ρ²、σ₅ | 定量断言 |
| `_verify_part5.py` | 基础 | 高斯波束、SSF 能量守恒（<1e-10）、PML、AOA=20° | 精确数值 |
| `_verify_part6.py` | 基础 | P-D 谱峰值、主模式映射、完整 pipeline | 端到端 |
| `_verify_radar.py` | Phase 1 | tau_2way_ms=2×tau_ms（误差<0.001 ms）、雷达 Pr<<Friis Pr（~126 dB）、phi_deg 正确 | 全部通过 |
| `_verify_OX.py` | Phase 2 | fH->0 收敛各向同性（差<1e-4）、12 O/X vs 6 iso 模式、tau 差<0.5 ms | 全部通过 |
| `_verify_absorption.py` | Phase 3 | A=0(chi=90)、f 依赖(1/f²)、beta 依赖(1/sinbeta)、雷达 2x 吸收（比值精确=2.00） | 全部通过 |
| `_verify_spreadf.py` | Phase 3 | Cs=0 无扰动、扰动集中在 h_screen、FFT 斜率~-(p+1)=-4（实测-3.21）、单调性 | 全部通过 |

---

## 9. 常见修改方法

### 改变工作频率

```python
# config.py
FREQ_MHZ = 15.0
# PE 的 dz_m 自动随 lambda 调整（_LAM_M / 4）
# 同步更新 GEOMAG: GEOMAG.update(get_geomag(IRI_LAT, IRI_LON, dt=IRI_DT))
```

### 改变 TX-RX 链路

```python
# config.py
TX_LAT   = 25.0
TX_LON   = 105.0
RX_RANGE = 800.0
RX_POS   = (800.0, 0.0)
IRI_LAT  = 27.0       # 路径中点纬度（约 TX + RX 纬度均值）
BG_X_MAX = 950.0      # 应大于 RX_RANGE + 200 km 余量
BG_X = np.arange(BG_X_MIN, BG_X_MAX + BG_DX, BG_DX)
```

### 切换日期/季节/昼夜

```python
# config.py
from datetime import datetime
IRI_DT = datetime(2021, 12, 1, 2, 0)   # 冬季夜间
# 夜间 D 层消失 -> 可开启 ABSORPTION + chi_deg=90 验证无吸收
```

### TID 参数敏感性分析

```python
import config as cfg
from models.hybrid_model import HybridPropagationModel

p2p_fast = {**cfg.P2P, 'n_init': 6, 'max_iter': 100}
results = []
for amp in [0.05, 0.10, 0.20, 0.30]:
    model = HybridPropagationModel(
        {'iri_params': {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
         'tid_params': {**cfg.TID, 'enable': True, 'amplitude': amp}},
        {'freq_MHz': cfg.FREQ_MHZ, 'Pt_W': cfg.PT_W, 'Gt': cfg.GT, 'Gr': cfg.GR},
    )
    modes, tau_ax, pd_W, main = model.compute(p2p_params=p2p_fast)
    results.append({'amp': amp, 'n_modes': len(modes), 'main': main})
```

### 调整 Spread-F 强度

```python
# Cs 越大，Ne 扰动越强，P-D 谱越展宽
spread_f_params = {**cfg.SPREAD_F, 'enable': True, 'Cs': 1e-2, 'p': 3.5}
```

### 比对日间与夜间 D 层吸收

```python
# 日间（chi=60°）
m_day = HybridPropagationModel(iono, rp,
    absorption_params={'enable': True, 'A0': 500.0, 'chi_deg': 60.0})

# 夜间（chi=90° -> A=0，无吸收）
m_night = HybridPropagationModel(iono, rp,
    absorption_params={'enable': True, 'A0': 500.0, 'chi_deg': 90.0})
```

### 实测 P-D 谱比对

```bash
# 准备测量数据 CSV（tau_ms, power_dBW 列）
# 先运行模型场景（生成 output/modes_baseline.csv）
conda run -n pytorch_cpu python -c "from main import run_baseline; run_baseline()"

# 执行比对
conda run -n pytorch_cpu python analysis/compare_pd.py \
    --measured  data/pd_measured_20240115.csv \
    --scenario  baseline \
    --tau_tol   0.5
```

### 更新 GEOMAG 至新路径中点

```python
from utils import get_geomag
from datetime import datetime

new_geomag = get_geomag(lat=32.5, lon=120.0, alt_km=300.0,
                         dt=datetime(2020, 1, 1))
# {'fH_MHz': 1.197, 'dip_deg': 48.7, 'decl_deg': -5.5}
# 将结果填入 config.py GEOMAG 字典
```

---

## 10. 性能参考

| 模块 | 默认参数 | 单次耗时（单线程估算） | 快速模式 |
|------|---------|----------------------|---------|
| M1 密度场 | IRI + TID | ~1 s | — |
| M2 P2P（各向同性） | n_init=18, max_iter=500 | ~30–60 s | n_init=6, max_iter=100: ~5 s |
| M2 P2P（O/X，2倍） | n_init=18, max_iter=500 | ~60–120 s | n_init=6: ~10 s |
| M3 Es | — | <0.1 s | — |
| M3b Spread-F | Nx=271 | <0.5 s | — |
| M4 PE/SSF | 100 km 域, dz=7.5 m | ~10–100 s | — |
| 完整 pipeline | M1+M2+M3+M4+M5-M7 | ~1–5 min | 仅开启 M2, n_init=6: ~10 s |

**并行加速**（P2P 为计算瓶颈）：
```python
# 调用脚本必须有 if __name__ == '__main__': 保护（Windows spawn 要求）
p2p = {**cfg.P2P, 'n_workers': 0}   # 0=自动检测 CPU 数
```

---

## 11. 参考文献

| # | 文献 | 对应模块 |
|---|------|---------|
| [1] | Hao, Y. et al. (2017). Sporadic-E layer scintillation. *Radio Science*. | M3 Es 三段式（fr=0.25, fs=0.10） |
| [2] | Koval, A. et al. (2018). TID modeling. | M1 TID（Hooke 1968 扰动公式） |
| [3] | Carrano, C. (2020). Wide-angle PE propagation. | M4 PE/SSF（宽角分步傅里叶） |
| [4] | Nosikov, I. et al. (2020). Variational P2P solver. *Radio Science*. | M2b 广义力变分 P2P 求解器 |
| [5] | Budden, K.G. (1961). *Radio Waves in the Ionosphere*. Cambridge. | M2a Appleton-Hartree O/X 公式 |
| [6] | Pederick, L.H. & Cervera, M.A. (2014). HF absorption. *Radio Science*, 49, 81–93. doi:10.1002/2013RS005274 | D 层吸收（A0=500 dB·MHz²） |
| [7] | Rino, C.L. (1979). Power law phase screen. *Radio Science*, 14(6), 1135–1145. | M3b SpreadFModel 幂律谱 |
| [8] | Ding, F. et al. (2021). HF through F-region irregularities. *Radio Science*. doi:10.1029/2020RS007239 | M3b 多层相位屏现代实现 |
| [9] | Jiang, C. et al. (2020). Plasma bubble modeling. | M1 等离子体泡物理参数 |
