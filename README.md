# HF Shortwave Hybrid Propagation Model

研究项目 李慧敏-22：复杂电离层环境短波反射/散射混合传播机理及建模研究

**测试链路**：TX @ 30°N 120°E，目标距离 1169 km，f = 10 MHz  
**代码仓库**：`https://github.com/zzmnl0/hf-propagation-model`  
**最新提交**：Phase 4 全部完成（Coleman 1997/1998 流管射线追踪 OTH 雷达后向散射）

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

# 单独运行流管雷达场景
conda run -n pytorch_cpu python -c "import sys,os;sys.path.insert(0,'.');from main import run_tube_radar;run_tube_radar()"

# 验证 Phase 1-4
conda run -n pytorch_cpu python tests/_verify_part6.py
conda run -n pytorch_cpu python tests/_verify_radar.py
conda run -n pytorch_cpu python tests/_verify_absorption.py
conda run -n pytorch_cpu python tests/_verify_spreadf.py
conda run -n pytorch_cpu python tests/_verify_tube.py

# 模式特征分析（需先运行 run_tube_radar()）
conda run -n pytorch_cpu python analysis/mode_analyzer.py --csv output/modes_tube_radar.csv

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
├── main.py                      # 场景入口（9个场景函数 + __main__）
│
├── models/
│   ├── ionosphere_model.py      # M1: IRI + TID / Es / 等离子体泡 / Spread-F
│   ├── ray_tracer.py            # M2a: Haselgrove RK4 + RefractiveIndexAH (O/X)
│   ├── point_to_point.py        # M2b: Nosikov 2020 变分 P2P 求解器
│   ├── tube_tracer.py           # M2c: Coleman 1997 流管射线追踪（Phase 4 新增）
│   ├── es_model.py              # M3: Hao 2017 Es 三段式模型
│   ├── pe_propagator.py         # M4: Carrano 2020 PE/SSF 等离子体泡散射
│   ├── spread_f_model.py        # M3b: Rino 1979 幂律相位屏（扩展 F）
│   └── hybrid_model.py          # M5-M7: 完整 pipeline（P-D 谱 + 主模式识别）
│
├── analysis/
│   ├── compare_pd.py            # 实测 vs 模型 P-D 谱比对框架（CLI 工具）
│   └── mode_analyzer.py         # 模式特征分析工具（Phase 4 新增，CLI 工具）
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
│   ├── _verify_spreadf.py       # Phase 3: 扩展 F 相位屏验证
│   └── _verify_tube.py          # Phase 4: 流管追踪 5 项验证（新增）
│
└── output/                      # 生成 PNG / CSV（gitignored）
```

---

## 3. 整体数据流

两条并行路径，由 `tube_mode` 开关选择：

```
config.py  （所有参数集中于此）
    |
    v
IonosphereModel.build_Ne_field(x, z)
    |-- 1. IRI-2016 一维背景 profile     [iri2016]
    |-- 2. TID 扰动  (Hooke 1968)        [enable=True]
    |-- 3. Es 薄层   (Hao 2017)          [enable=True]
    |-- 4. 等离子体泡 (Gaussian 耗散)    [enable=True]
    |-- 5. Spread-F  (Rino 1979 相位屏)  [enable=True]
    |
    v  Ne_2d (Nx, Nz) [m^-3]
    |
    |══════════════════════════════════════════════════════════════════
    |   tube_mode=False（默认，P2P 变分路径）
    |══════════════════════════════════════════════════════════════════
    |
    v
RefractiveIndex（各向同性）或 RefractiveIndexAH（Appleton-Hartree O/X）
    |
    v
find_all_rays_p2p(tx, rx, n_model, freq, wave_mode='both'|None)
    |-- wave_mode=None  -> 各向同性，wave_mode='iso'
    |-- wave_mode='both' -> 分别用 O/X 跑变分，label 追加 _O/_X
    |
    v  rays: list[dict]  (points, tau_ms, h_reflect_km, wave_mode, label, ...)
    |
    v
对每条 ray（Power pipeline）：
    |-- 通信模式: Pr = Pt*Gt*Gr / L_free（Friis）
    |-- 雷达模式: Pr = Pt*Gt*Gr*lambda^2*sigma / ((4pi)^3 * R^4)
    |-- D 层吸收修正（可选）: A_dB = A0*cos(chi)^0.75 / ((f+fH_L)^2*sin(beta))
    |-- Es 修正（可选）: Pr_W 修正 + delta_tau_ms 累加
    |-- PE/SSF 泡散射（可选）: Pr_W *= P_out/P_in + delta_tau_ms
    |
    |══════════════════════════════════════════════════════════════════
    |   tube_mode=True（流管路径，需同时 radar_mode=True）
    |══════════════════════════════════════════════════════════════════
    |
    v
TubeRayTracer.compute(tx, x_tgt, Pt, Gt, Gr, sigma0)
    |-- A. shoot_fan(): 扇形扫射 beta_min~beta_max，步长 delta_beta_deg
    |-- B. compute_tubes(): 相邻射线对形成流管
    |       F_focus = R_eff * d_beta / d_x_land   （聚焦因子）
    |       A_tube  = delta_x * L_cross_km         （地面投影面积）
    |-- C. newton_refine(): Newton 迭代精化落点
    |       F(beta) = x_land(beta) - x_tgt = 0
    |       dF/dbeta 用中心差分（d_beta=0.01 deg）
    |-- D. backscatter_power_W() + pulse_correction()
    |       Pr = Pt*Gt*Gr*lambda^2*sigma0*A_tube / ((4pi)^3 * P_one^4)
    |       f_pulse = T_pulse/tau_spread  if tau_spread > T_pulse else 1
    |-- E. _dedup_modes(): 同标签 + |dtau| < 0.05 ms -> 保留最大功率
    |-- 同样支持 D 层吸收 + Es 修正（两路统一处理）
    |
    |══════════════════════════════════════════════════════════════════
    |   两路汇合
    |══════════════════════════════════════════════════════════════════
    |
    v  mode_results: list[dict]
    |
    v
build_pd_spectrum(mode_results, method='gaussian'|'tube')
    |-- gaussian: Pr_i * exp(-(tau-tau_i)^2/(2*sigma_i^2))
    |-- tube:     矩形窗（宽度=tau_spread_ms），物理反映管内时延扩展
    |
    v
identify_main_mode(tau_axis, pd_W, mode_results)
    -> (main_mode, ranked_modes, mode_summary)
    mode_summary: {n_modes, main_label, main_tau_ms, main_Pr_dBW,
                   main_F_focus, tau_spread_main, mode_pairs}
    |
    v
model.compute() 返回: (mode_results, tau_axis, pd_W, main_mode)
model.mode_summary                    <- 额外实例属性（Phase 4 新增）
```

---

## 4. 模块详解

### M1 — IonosphereModel

**文件**：`models/ionosphere_model.py`

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
n^2_{O,X} = 1 - X / (1 - YT^2/[2(1-X)] -/+ sqrt(YT^4/[4(1-X)^2] + YL^2))

X  = fp^2/f^2          YT = Y*sin(alpha)   alpha = pi/2 - dip
Y  = fH/f              YL = Y*cos(alpha)
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

```python
rays = find_all_rays_p2p(
    tx_km, rx_km, n_model, freq_MHz,
    p2p_params = P2P,
    wave_mode  = 'both',   # None=各向同性 | 'O' | 'X' | 'both'
    geomag     = GEOMAG,
)
```

**算法**（Nosikov 2020）：

- 光程泛函 `S = sum n(mid_i) * |seg_i|`（离散 Fermat 积分）
- 梯度下降 `pts -= alpha * grad_perp` + 弹簧力平滑
- 去重：h 差 < `clust_h_km` 且 tau 差 < `clust_tau_ms`

**模式标签规则**：

| h_reflect | tau_ms | 基础标签 | O/X 开启后 |
|-----------|--------|----------|-----------|
| < 140 km | — | `Es` | `Es_O` / `Es_X` |
| 140–200 km | — | `E` | `E_O` / `E_X` |
| 200–300 km | < 5 ms | `1F_low` | `1F_low_O` / `1F_low_X` |
| 200–300 km | >= 5 ms | `1F_high` | `1F_high_O` / ... |
| >= 300 km | — | `2F` | `2F_O` / `2F_X` |

各向同性时 wave_mode = `'iso'`，标签无后缀。

**辅助函数**：

| 函数 | 说明 |
|------|------|
| `classify_mode(ray_dict)` | 返回含 O/X 后缀的标签 |
| `extract_es_params(points, h_Es_km)` | 射线与 Es 层的交叉参数 |
| `extract_bubble_entry(points, z_bot_km)` | 射线进入等离子体泡的参数 |

---

### M2c — TubeRayTracer（Phase 4 新增）

**文件**：`models/tube_tracer.py`  **类**：`TubeRayTracer`

OTH 雷达后向散射物理模型（Coleman 1997/1998）。用相邻射线对（流管）替代简单雷达方程，计算包含聚焦/散焦效应的后向散射功率。

```python
from models.tube_tracer import TubeRayTracer

tracer = TubeRayTracer(
    n_model     = nm,           # RefractiveIndex 实例
    freq_MHz    = 10.0,
    tube_params = TUBE_TRACER,  # 可选，默认从 config 读取
    rt_params   = RT,
)

sigma0 = 10.0 ** (RADAR['sigma0_ground_dB'] / 10.0)   # 线性地面 RCS

tube_modes = tracer.compute(
    tx_pos  = (0.0, 0.0),   # 发射机位置 [km]
    x_tgt   = 1169.0,        # 目标水平距离 [km]
    Pt_W    = 1000.0,
    Gt      = 1.0,
    Gr      = 1.0,
    sigma0  = sigma0,
    newton  = True,          # 开启 Newton 落点精化
)
```

**核心公式**：

```
# Coleman (1997, Radio Science 32(1), Eq.15-17)
F_focus    = R_eff * d_beta / d_x_land     [无量纲，>1=聚焦，<1=散焦]
A_tube     = delta_x * L_cross_km          [km^2，地面投影面积]

# 后向散射功率（与 radar_equation_W 公式等价）
Pr = Pt*Gt*Gr * lambda^2 * sigma0 * A_tube
     / ((4*pi)^3 * P_one^4)

# 脉冲修正（管内时延扩展超过脉宽时）
f_pulse = T_pulse / tau_spread   if tau_spread > T_pulse else 1.0
```

**关键参数说明**：

- `R_eff` = 单程群路径（group_path_km），约 1200 km（F 层典型）
- `delta_x` = 相邻射线落点距离差 [km]
- `L_cross_km` = 方位向波束足迹（2D->3D 面积扩展）
- Newton 精化：F(beta) = x_land(beta) - x_tgt = 0，中心差分 Jacobian

**方法一览**：

| 方法 | 说明 |
|------|------|
| `shoot_fan(tx_pos)` | 扇形扫射，返回带 x_land 的射线列表 |
| `compute_tubes(fan_rays, x_tgt)` | 识别落点近 x_tgt 的相邻射线对，计算流管几何 |
| `backscatter_power_W(Pt, Gt, Gr, P_one_km, sigma0, A_tube_km2)` | Coleman 功率公式 |
| `pulse_correction(tau_spread_ms)` | 脉冲展宽修正因子 |
| `newton_refine(tx_pos, beta0, x_tgt)` | Newton 迭代精化中心射线，返回 (ray, beta_final) |
| `compute(...)` | 完整流管 pipeline，返回 mode_results 列表 |

**通过 HybridPropagationModel 调用**（推荐，自动集成 D 层吸收 + Es 修正）：

```python
model = HybridPropagationModel(
    iono_params  = _make_iono_params(),
    radar_params = {
        'freq_MHz'        : 10.0,
        'Pt_W'            : 1000.0,
        'Gt'              : 1.0,
        'Gr'              : 1.0,
        'sigma_rcs_m2'    : 5.0,
        'sigma0_ground_dB': -20.0,   # 地面归一化 RCS [dBsm]
    },
    radar_mode = True,
    tube_mode  = True,    # 开启流管模式（同时需要 radar_mode=True）
)
modes, tau_ax, pd, main = model.compute(TX_POS, RX_POS)
summary = model.mode_summary   # 见 §6 输出格式
```

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

```python
sfm = SpreadFModel(Cs=1e-3, p=3.0, h_screen_km=300.0, L0_km=50.0, seed=42)
Ne_2d_new = sfm.apply(Ne_2d, x_km, z_km)
```

功率谱 `S(k) ~ (k²+k₀²)^(-(p+1)/2)`，垂直 Gaussian 包络（半宽 20 km）。

---

### M4 — PE/SSF 传播器

**文件**：`models/pe_propagator.py`  **类**：`PEPropagator`

SSF 单步（Carrano 2020 宽角 PE）：
- Step A（折射，空域）：`u_half = u * exp(i k0 (n-1) dx)`
- Step B（衍射，谱域）：FFT -> kx_eff -> IFFT

地球展平修正（`PE['earth_flat']=True`）：`n_eff = n * (1 + z/R_E)`

---

### M5-M7 — HybridPropagationModel

**文件**：`models/hybrid_model.py`

```python
model = HybridPropagationModel(
    iono_params       = {...},        # 见下方结构
    radar_params      = {...},        # freq_MHz, Pt_W, Gt, Gr [, sigma0_ground_dB]
    radar_mode        = False,        # True = 雷达方程 + tau_2way_ms
    geomag_params     = None,         # {**GEOMAG, 'enable_OX': True} = O/X 分裂
    absorption_params = None,         # {enable:True, A0:500, ...} = D 层吸收
    tube_mode         = False,        # True = 流管追踪（需 radar_mode=True）
)
modes, tau_ax, pd_W, main = model.compute(tx_km, rx_km, t=0.0, p2p_params=P2P)
summary = model.mode_summary         # Phase 4 新增实例属性
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
# tube_mode=False，radar_mode=False（通信）
Pr_W = Pt * Gt * Gr * 10^(-L_free_dB/10)

# tube_mode=False，radar_mode=True（雷达方程）
Pr_W = Pt*Gt*Gr*lambda^2*sigma / ((4pi)^3 * R^4)
tau_2way_ms = 2 * tau_ms

# tube_mode=True，radar_mode=True（Coleman 流管）
Pr_W = Pt*Gt*Gr*lambda^2*sigma0*A_tube / ((4pi)^3 * P_one^4)
tau_2way_ms = 2 * tau_ms
# 附加：f_pulse = T_pulse/tau_spread（脉冲展宽修正）
# 附加：F_focus = R_eff*d_beta/d_x_land（聚焦因子）
```

**P-D 谱构建**（Phase 4 新增 tube 模式）：

```python
# method='gaussian'（默认，P2P 路径）
# sigma_i = max(delta_tau_ms_i, tau_res_ms)

# method='tube'（流管路径）
# 矩形窗宽度 = max(tau_spread_ms, delta_tau_ms, tau_res_ms)
# 物理意义：管内时延扩展均匀分布
```

**模块级函数**：

| 函数 | 签名 | 说明 |
|------|------|------|
| `build_pd_spectrum` | `(mode_results, tau_axis, tau_res_ms, method)` | 构建 P-D 谱 |
| `identify_main_mode` | `(tau_axis, pd_W, mode_results)` | 峰值检测，返回 3 元组 |
| `_find_ox_pairs` | `(mode_results)` | O/X 配对匹配 |

---

### 辅助工具

**`utils.py`** 中的主要函数：

| 函数 | 说明 |
|------|------|
| `free_space_loss_dB(D_km, freq_MHz)` | Friis 自由空间损耗 [dB] |
| `radar_equation_W(Pt, Gt, Gr, freq, gp_km, sigma)` | 单基雷达接收功率 [W] |
| `d_layer_absorption_dB(freq, beta, chi, A0, fH_L)` | D 层非偏吸收 [dB] |
| `to_dBW(Pr_W)` | 线性功率 -> dBW |
| `get_geomag(lat, lon, alt_km, dt)` | ppigrf IGRF-14 地磁参数 |
| `haversine_km(lat1,lon1,lat2,lon2)` | 大圆距离 [km] |
| `bearing_deg(lat1,lon1,lat2,lon2)` | 初始方位角 [deg] |

**`analysis/compare_pd.py`** 实测比对工具（CLI）：

```bash
python analysis/compare_pd.py \
    --measured  data/pd_measured.csv \
    --scenario  baseline \
    --tau_tol   0.5 \
    --output    output/compare_baseline.png
```

**`analysis/mode_analyzer.py`** 模式特征分析工具（Phase 4 新增，CLI）：

```bash
# 分析模式特征并打印报告
python analysis/mode_analyzer.py --csv output/modes_tube_radar.csv

# 也可作为库调用
from analysis.mode_analyzer import (analyze_ox_pairs, analyze_mode_features,
                                     build_mode_summary, print_mode_report)

features  = analyze_mode_features(mode_results)   # 时延/功率差矩阵
ox_pairs  = analyze_ox_pairs(mode_results)         # O/X 配对（P2P 路径）
summary   = build_mode_summary(mode_results)        # 与 model.mode_summary 等价
print_mode_report(mode_results, summary)            # 格式化输出
```

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
| `sigma_rcs_m2` | 5.0 | 目标 RCS [m²]（10 MHz 飞机谐振区） |
| `two_way` | True | 输出双程时延 tau_2way_ms |
| `sigma0_ground_dB` | -20.0 | 地面归一化 RCS [dBsm]（海面 -20，陆地 -25） |
| `LINK_BEARING_DEG` | 0.0 | TX → 目标方位角 [°]（正北=0） |

### 流管追踪参数（`config.TUBE_TRACER`，Phase 4 新增）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `delta_beta_deg` | 0.5 | 扇形扫射仰角步长 [°] |
| `n_tube_rays` | 80 | 最大扫射射线数 |
| `x_tgt_tol_km` | 80.0 | 目标落点搜索窗口半宽 [km] |
| `T_pulse_ms` | 0.5 | 雷达脉冲宽度 [ms]（脉冲修正参考） |
| `newton_tol_km` | 1.0 | Newton 精化收敛阈值 [km] |
| `newton_max_iter` | 5 | Newton 最大迭代次数 |
| `L_cross_km` | 100.0 | 方位向波束足迹 [km]（2D->3D 面积扩展） |

### 地磁参数（`config.GEOMAG`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `fH_MHz` | 1.197 | 回旋频率 [MHz]（B=42766 nT，ppigrf IGRF-14） |
| `dip_deg` | 48.7 | 地磁倾角 [°] |
| `decl_deg` | -5.5 | 地磁偏角 [°]（负=偏西） |
| `enable_OX` | False | O/X 磁离子分裂总开关 |

更新方式：`from utils import get_geomag; GEOMAG.update(get_geomag(lat, lon, dt=dt))`

### D 层吸收参数（`config.ABSORPTION`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | D 层吸收修正开关 |
| `A0` | 500.0 | 吸收系数 [dB·MHz²]（Pederick & Cervera 2014） |
| `chi_deg` | 60.0 | 太阳天顶角 [°]（60=白天典型，90=夜间=0 dB） |

### TID 参数（`config.TID`）

| 键 | 默认值 | 说明 |
|-----|--------|------|
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
| `store_history` | False | 存储所有 x 截面（调试，~GB RAM） |

---

## 6. 输出格式

### mode_results 字典键

| 键 | 单位 | 说明 | 来源 |
|----|------|------|------|
| `label` | — | 模式标签（见分类规则） | 两路 |
| `tau_ms` | ms | 单程群时延 | 两路 |
| `tau_2way_ms` | ms | 双程时延 | 两路（雷达/流管模式） |
| `delta_tau_ms` | ms | Es/PE 引起的时延扩展 | P2P 路径 |
| `tau_spread_ms` | ms | 管内时延扩展（相邻射线时延差） | 流管路径 |
| `Pr_W` | W | 接收功率（线性） | 两路 |
| `Pr_dBW` | dBW | 接收功率（对数） | 两路 |
| `h_reflect_km` | km | 反射高度（射线最高点） | 两路 |
| `group_path_km` | km | 单程群路径长度 | 两路 |
| `beta_deg` | deg | 出发仰角 | 两路 |
| `phi_deg` | deg | 链路方位角（来自 LINK_BEARING_DEG） | 两路 |
| `wave_mode` | — | `'O'` / `'X'` / `'iso'` | 两路 |
| `F_focus` | — | 聚焦因子（>1=聚焦，<1=散焦） | 流管路径 |
| `A_tube_km2` | km² | 流管地面投影面积 | 流管路径 |
| `points` | km | 控制点数组 (n_ctrl+2, 2) | P2P 路径 |

### mode_summary 字典（Phase 4 新增，`model.mode_summary`）

| 键 | 类型 | 说明 |
|----|------|------|
| `n_modes` | int | 传播模式总数 |
| `main_label` | str | 主模式标签 |
| `main_tau_ms` | float | 主模式单程时延 [ms] |
| `main_Pr_dBW` | float | 主模式功率 [dBW] |
| `main_F_focus` | float | 主模式聚焦因子（P2P 路径时为 1.0） |
| `tau_spread_main` | float | 主模式时延扩展 [ms] |
| `mode_pairs` | list[dict] | O/X 配对列表（各向同性时为空列表） |

`mode_pairs` 每个元素的键：`label_O, label_X, tau_O, tau_X, delta_tau_OX_ms, Pr_O_dBW, Pr_X_dBW, delta_Pr_OX_dB`

### identify_main_mode 返回值（Phase 4 更新为 3 元组）

```python
main_mode, ranked_modes, mode_summary = identify_main_mode(tau_axis, pd_W, mode_results)
# main_mode    : dict | None   P-D 谱主峰对应的模式
# ranked_modes : list[dict]    按 P-D 峰值功率降序排列的模式列表
# mode_summary : dict          见上表
```

### compute() 返回值

```python
modes, tau_ax, pd_W, main = model.compute(tx_km, rx_km)
# modes   : list[dict]         各传播模式（含所有键）
# tau_ax  : np.ndarray [ms]    P-D 谱横轴
# pd_W    : np.ndarray [W]     P-D 谱功率
# main    : dict | None        主模式
# model.mode_summary           额外属性（随 compute() 更新）
```

### CSV 输出列

`save_modes_csv(scenario, modes, tau_ax, pd_W, freq_MHz, out_dir)` 输出固定列：

```
scenario, label, freq_MHz, tau_ms, tau_2way_ms, delta_tau_ms,
Pr_W, Pr_dBW, pd_power_W, h_reflect_km, group_path_km, beta_deg, phi_deg
```

> 流管新字段（F_focus, A_tube_km2, tau_spread_ms）不写入 CSV，由 `model.mode_summary` 和 `print_mode_report()` 输出。

---

## 7. 场景运行示例

### 当前 9 个场景（main.py）

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
| `run_tube_radar()` | IRI | Coleman 流管后向散射 | `pd_tube_radar`（Phase 4 新增） |

### 代码示例

```python
import config as cfg
from models.hybrid_model import HybridPropagationModel

# ── 最简运行（IRI 仅背景，通信模式）──────────────────────────────────────────
model = HybridPropagationModel(
    iono_params  = {'iri_params': {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON}},
    radar_params = {'freq_MHz': 10.0, 'Pt_W': 1000.0, 'Gt': 1.0, 'Gr': 1.0},
)
modes, tau_ax, pd_W, main = model.compute()

# ── OTH 雷达方程模式 ──────────────────────────────────────────────────────────
model = HybridPropagationModel(
    iono_params  = {'iri_params': {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON}},
    radar_params = {'freq_MHz': 10.0, 'Pt_W': 1e6, 'Gt': 1.0, 'Gr': 1.0,
                    'sigma_rcs_m2': cfg.RADAR['sigma_rcs_m2']},
    radar_mode   = True,
)

# ── OTH 雷达流管模式（Coleman 1997，物理聚焦）────────────────────────────────
model = HybridPropagationModel(
    iono_params  = {'iri_params': {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON}},
    radar_params = {'freq_MHz': 10.0, 'Pt_W': 1000.0, 'Gt': 1.0, 'Gr': 1.0,
                    'sigma0_ground_dB': cfg.RADAR['sigma0_ground_dB']},
    radar_mode   = True,
    tube_mode    = True,
)
modes, tau_ax, pd_W, main = model.compute(cfg.TX_POS, cfg.RX_POS)
summary = model.mode_summary
# summary['main_F_focus'] -> 聚焦因子（测试链路典型值 ~0.45，散焦）
# summary['tau_spread_main'] -> 管内时延扩展 [ms]

# ── O/X 磁离子分裂 ────────────────────────────────────────────────────────────
model = HybridPropagationModel(
    iono_params   = {'iri_params': {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
                     'tid_params': {**cfg.TID, 'enable': True}},
    radar_params  = {'freq_MHz': 10.0, 'Pt_W': 1000.0, 'Gt': 1.0, 'Gr': 1.0},
    geomag_params = {**cfg.GEOMAG, 'enable_OX': True},
)

# ── 开启 D 层吸收修正 ─────────────────────────────────────────────────────────
model = HybridPropagationModel(
    iono_params       = {...},
    radar_params      = {...},
    absorption_params = {**cfg.ABSORPTION, 'enable': True},
)

# ── 快速测试（减少 P2P 迭代）─────────────────────────────────────────────────
p2p_fast = {**cfg.P2P, 'n_init': 6, 'max_iter': 100}
modes, *_ = model.compute(p2p_params=p2p_fast)

# ── P2P 并行加速（脚本必须有 __main__ 保护）──────────────────────────────────
if __name__ == '__main__':
    modes, *_ = model.compute(p2p_params={**cfg.P2P, 'n_workers': 0})
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
| `_verify_radar.py` | Phase 1 | tau_2way=2×tau（误差<0.001 ms）、Pr_radar<<Pr_friis（~126 dB）、phi_deg 正确 | 全部通过 |
| `_verify_OX.py` | Phase 2 | fH->0 收敛各向同性（差<1e-4）、O/X 模式数=2×iso、tau 差<0.5 ms | 全部通过 |
| `_verify_absorption.py` | Phase 3 | A=0(chi=90)、f 依赖(1/f²)、beta 依赖(1/sin)、雷达 2x 吸收（比值精确=2.00） | 全部通过 |
| `_verify_spreadf.py` | Phase 3 | Cs=0 无扰动、扰动集中 h_screen、FFT 斜率~-(p+1)=-4、单调性 | 全部通过 |
| `_verify_tube.py` | Phase 4 | tau 与 P2P 差<0.05 ms、功率差 0.00 dB、F_focus 范围、TID 增大 tau_spread、mode_summary 完整 | 全部通过 |

**Phase 4 验证数值（IRI baseline，1169 km，10 MHz）**：

| Check | 内容 | 结果 |
|-------|------|------|
| 1 | 1F_low: tube=4.214 ms vs P2P=4.202 ms | diff=0.013 ms < 0.05 ms |
| 2 | 功率一致性（同 sigma） | diff=0.00 dB |
| 3 | F_focus=0.454 | in (0.1, 10.0) |
| 4 | TID: tau_spread 0.078->0.097 ms | 单调增加 |
| 5 | mode_summary: n_modes=1, main=1F_low, F_focus=0.45 | 完整 |

---

## 9. 常见修改方法

### 切换流管 / 雷达方程模式

```python
# 雷达方程模式（近似，快）
model = HybridPropagationModel(iono_params, radar_params, radar_mode=True)

# 流管模式（物理，含聚焦因子，慢约 3-5x）
model = HybridPropagationModel(iono_params, radar_params,
                               radar_mode=True, tube_mode=True)
```

### 调整地面 RCS（影响流管模式功率）

```python
# config.py
RADAR['sigma0_ground_dB'] = -25.0   # 陆地（默认 -20.0=海面）

# 也可在 radar_params 中直接覆盖
radar_params = {'freq_MHz': 10.0, ..., 'sigma0_ground_dB': -30.0}
```

### 调整流管密度（影响模式分辨率）

```python
# config.py
TUBE_TRACER['delta_beta_deg'] = 0.2   # 更细分辨率（慢）
TUBE_TRACER['x_tgt_tol_km']  = 50.0  # 收紧搜索窗口
TUBE_TRACER['newton_max_iter'] = 8    # 提高精化精度
```

### 改变工作频率

```python
# config.py
FREQ_MHZ = 15.0
# PE 的 dz_m 自动随 lambda 调整（_LAM_M / 4）
# 更新地磁参数: GEOMAG.update(get_geomag(IRI_LAT, IRI_LON, dt=IRI_DT))
```

### 改变 TX-RX 链路

```python
# config.py
TX_LAT   = 25.0
TX_LON   = 105.0
RX_RANGE = 800.0
RX_POS   = (800.0, 0.0)
RADAR['target_range_km'] = 800.0    # 流管模式目标距离同步修改
IRI_LAT  = 27.0
BG_X_MAX = 950.0    # 大于 RX_RANGE + 200 km
BG_X = np.arange(BG_X_MIN, BG_X_MAX + BG_DX, BG_DX)
```

### TID 参数敏感性分析

```python
import config as cfg
from models.hybrid_model import HybridPropagationModel

p2p_fast = {**cfg.P2P, 'n_init': 6, 'max_iter': 100}
for amp in [0.05, 0.10, 0.20, 0.30]:
    model = HybridPropagationModel(
        {'iri_params': {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
         'tid_params': {**cfg.TID, 'enable': True, 'amplitude': amp}},
        {'freq_MHz': cfg.FREQ_MHZ, 'Pt_W': cfg.PT_W, 'Gt': cfg.GT, 'Gr': cfg.GR},
    )
    modes, _, _, _ = model.compute(p2p_params=p2p_fast)
    print("amp={:.2f}: n_modes={}".format(amp, len(modes)))
```

### 流管模式 TID 影响分析

```python
# 流管模式下 TID 通过 tau_spread_ms 量化管内时延扩展
from models.tube_tracer import TubeRayTracer
from models.ionosphere_model import IonosphereModel
from models.ray_tracer import RefractiveIndex
import config as cfg

sigma0 = 10.0 ** (cfg.RADAR['sigma0_ground_dB'] / 10.0)

for enable in [False, True]:
    iono = IonosphereModel(tid_params={**cfg.TID, 'enable': enable})
    Ne, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    nm = RefractiveIndex(Ne, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ)
    tracer = TubeRayTracer(nm, cfg.FREQ_MHZ)
    modes = tracer.compute(cfg.TX_POS, cfg.RADAR['target_range_km'],
                           cfg.PT_W, cfg.GT, cfg.GR, sigma0)
    for m in modes:
        print("TID={}: {} tau_spread={:.4f}ms".format(enable, m['label'], m['tau_spread_ms']))
```

### 比对日间与夜间 D 层吸收

```python
m_day = HybridPropagationModel(iono, rp,
    absorption_params={'enable': True, 'A0': 500.0, 'chi_deg': 60.0})

m_night = HybridPropagationModel(iono, rp,
    absorption_params={'enable': True, 'A0': 500.0, 'chi_deg': 90.0})
```

### 实测 P-D 谱比对

```bash
conda run -n pytorch_cpu python -c "import sys,os;sys.path.insert(0,'.');from main import run_tube_radar;run_tube_radar()"

conda run -n pytorch_cpu python analysis/compare_pd.py \
    --measured  data/pd_measured.csv \
    --scenario  tube_radar \
    --tau_tol   0.5
```

---

## 10. 性能参考

| 模块 / 场景 | 配置 | 耗时估算（单线程） | 快速模式 |
|-------------|------|-------------------|---------|
| M1 密度场 | IRI + TID | ~1 s | — |
| M2b P2P（各向同性） | n_init=18, max_iter=500 | ~30–60 s | n_init=6, iter=100: ~5 s |
| M2b P2P（O/X，2x） | n_init=18, max_iter=500 | ~60–120 s | n_init=6: ~10 s |
| M2c 流管 | 80 射线，Newton=True | ~10–30 s | Newton=False: ~5 s |
| M3 Es | — | <0.1 s | — |
| M3b Spread-F | Nx=271 | <0.5 s | — |
| M4 PE/SSF | 100 km 域, dz=7.5 m | ~10–100 s | — |
| run_tube_radar() | IRI + 流管 + Newton | ~15–40 s | tube_params n_tube_rays=40: ~8 s |
| 完整 pipeline (P2P) | M1+M2b+M3+M4+M5-M7 | ~1–5 min | n_init=6: ~15 s |

**并行加速**（P2P 为计算瓶颈，流管单线程）：
```python
if __name__ == '__main__':
    p2p = {**cfg.P2P, 'n_workers': 0}   # Windows spawn 要求 __main__ 保护
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
| [6] | Pederick, L.H. & Cervera, M.A. (2014). HF absorption. *Radio Science*, 49, 81–93. | D 层吸收（A0=500 dB·MHz²） |
| [7] | Rino, C.L. (1979). Power law phase screen. *Radio Science*, 14(6), 1135–1145. | M3b SpreadFModel 幂律谱 |
| [8] | Ding, F. et al. (2021). HF through F-region irregularities. *Radio Science*. doi:10.1029/2020RS007239 | M3b 多层相位屏现代实现 |
| [9] | Jiang, C. et al. (2020). Plasma bubble modeling. | M1 等离子体泡物理参数 |
| [10] | Coleman, C.J. (1997). A ray tracing formulation and its application to some problems in over-the-horizon radar. *Radio Science*, 32(1), 45–60. | M2c TubeRayTracer 流管几何（F_focus 公式 Eq.15-17） |
| [11] | Coleman, C.J. (1998). A model of HF sky wave radar clutter. *Radio Science*, 33(4), 921–929. | M2c 后向散射功率与脉冲修正 |
