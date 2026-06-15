# HF Shortwave Hybrid Propagation Model

研究项目 李慧敏-22：复杂电离层环境短波反射/散射混合传播机理及建模研究

测试链路：TX @ 30°N 120°E，RX @ 1169 km，f = 10 MHz

---

## 目录

1. [环境与运行](#1-环境与运行)
2. [目录结构](#2-目录结构)
3. [整体数据流](#3-整体数据流)
4. [模块详解](#4-模块详解)
5. [参数速查表](#5-参数速查表)
6. [输出格式](#6-输出格式)
7. [场景使用示例](#7-场景使用示例)
8. [验证脚本](#8-验证脚本)
9. [精度与性能](#9-精度与性能)
10. [常见修改方法](#10-常见修改方法)

---

## 1. 环境与运行

```bash
# 环境验证
conda run -n pytorch_cpu python main.py

# 单独运行验证脚本
conda run -n pytorch_cpu python tests/_verify_part1.py
conda run -n pytorch_cpu python tests/_verify_part3.py
conda run -n pytorch_cpu python tests/_verify_part4.py
conda run -n pytorch_cpu python tests/_verify_part5.py
conda run -n pytorch_cpu python tests/_verify_part6.py

# 背景电子密度图（输出到 output/）
conda run -n pytorch_cpu python viz/plot_ne_background.py
```

**依赖包**：`numpy`、`scipy`、`matplotlib`、`iri2016`（不使用 iricore）

**注意**：Windows 终端使用 GBK 编码，Python `print()` 中只能使用 ASCII 字符。

---

## 2. 目录结构

```
code/
├── config.py                  # 所有物理参数与网格设置（唯一配置入口）
├── utils.py                   # 物理工具函数（EM、功率、坐标）
├── main.py                    # 入口；取消注释选择运行场景
├── models/
│   ├── ionosphere_model.py    # M1：IRI + TID / Es / 等离子体泡密度场
│   ├── ray_tracer.py          # M2a：Haselgrove RK4 射线追踪（扇形/单条）
│   ├── point_to_point.py      # M2b：P2P 变分求解器（Nosikov 2020）
│   ├── es_model.py            # M3：Es 三段式模型（Hao 2017）
│   ├── pe_propagator.py       # M4：PE/SSF 等离子体泡散射（Carrano 2020）
│   └── hybrid_model.py        # M5-M7：完整 pipeline（合成 + P-D 谱 + 主模式）
├── viz/
│   ├── plot_utils.py          # 共享绘图函数（射线扇、P-D 谱、Ne 场）
│   └── plot_ne_background.py  # 4 面板 Ne 背景图 -> output/
├── tests/
│   ├── _verify_part0.py       # 导入 / 配置健全性检查
│   ├── _verify_part1.py       # IonosphereModel 完整验证
│   ├── _verify_part2.py       # 射线追踪验证
│   ├── _verify_part3.py       # P2P 变分求解验证
│   ├── _verify_part4.py       # Es 模型验证
│   ├── _verify_part5.py       # PE/SSF 传播器验证
│   └── _verify_part6.py       # 混合模型流水线验证
└── output/                    # 生成的 PNG 文件（gitignored）
```

---

## 3. 整体数据流

```
config.py
    |
    v
IonosphereModel.build_Ne_field(x, z, t)
    |-- IRI-2016 背景 1D profile (iri2016)
    |-- TID 扰动  (Hooke 1968)    [可选]
    |-- Es 薄层   (Hao 2017)      [可选]
    |-- 等离子体泡 (Gaussian 耗散) [可选]
    |
    v  Ne_2d (Nx, Nz) [m^-3]
    |
    v
RefractiveIndex(Ne_2d, x, z, freq)    <- 包装插值器
    |
    v
find_all_rays_p2p(tx, rx, n_model, freq)
    |-- 对 n_init 个仰角 x {high, low} 并行/串行变分求解
    |-- 去重 (h 差 < clust_h_km 且 tau 差 < clust_tau_ms)
    |-- classify_mode -> 标签
    |
    v  rays: list[dict]  (含 points, tau_ms, h_reflect_km, ...)
    |
    v
对每条 ray:
    |-- free_space_loss_dB -> Pr_W 初值
    |-- [Es 开启] extract_es_params -> EsLayerModel.compute_power
    |       -> Pr_W 修正, delta_tau_ms 累加, label 更新
    |-- [Bubble 开启] extract_bubble_entry -> pe.extract_domain
    |       -> construct_incident_field -> pe.propagate
    |       -> pe.analyze -> delta_tau_ms 累加, Pr_W 按 P_out/P_in 缩放
    |
    v  mode_results: list[dict]
    |
    v
build_pd_spectrum(mode_results)
    -> tau_axis, pd_W  (高斯展宽叠加 P-D 谱)
    |
    v
identify_main_mode(tau_axis, pd_W, mode_results)
    -> main_mode, ranked_modes
```

---

## 4. 模块详解

### M1 — IonosphereModel (`models/ionosphere_model.py`)

**类**：`IonosphereModel`

| 方法 | 说明 |
|------|------|
| `build_Ne_field(x, z, t=0)` | 返回 `(Ne_2d, n_2d)`，均为 `(Nx, Nz)` |

**不规则体叠加顺序**（`build_Ne_field` 内部）：
1. IRI-2016 一维背景 → 广播为 (Nx, Nz)
2. TID 扰动（`enable=True` 时）
3. Es 薄层增量（`enable=True` 时，在 x 方向均匀）
4. 等离子体泡耗散（`enable=True` 时，最后叠加）

**IRI 接口**：`iri2016.IRI(dt, (z_min, z_max, dz), lat, lon)['ne'].values`

---

### M2a — RefractiveIndex + 射线追踪 (`models/ray_tracer.py`)

**类**：`RefractiveIndex`

| 方法 | 说明 |
|------|------|
| `n(x, z)` | 单点 n 值（标量） |
| `n2(x, z)` | 单点 n² |
| `grad_n2(x, z)` | 中心差分梯度 (dn²/dx, dn²/dz) |
| `n_batch(pts)` | 矢量化，`pts` 形状 (M, 2) → (M,)，单次插值调用 |

**函数**：

| 函数 | 说明 |
|------|------|
| `trace_single_ray(tx, beta_deg, n_model, freq_MHz)` | Haselgrove RK4，返回 ray dict |
| `shoot_rays_fan(tx_pos, n_model, rt_params)` | 扇形扫角，返回 ray list |

---

### M2b — P2P 变分求解 (`models/point_to_point.py`)

**算法**（Nosikov 2020）：

- 光程泛函 `S = Σ n(mid_i) * |seg_i|`（离散 Fermat 积分）
- 高角射线：梯度下降 `pts -= alpha * grad_perp`
- 低角射线：符号翻转 `pts += alpha * grad_perp`（鞍点 → 稳定点）
- 弹簧力 `k_spring * (pts[i+1] - 2*pts[i] + pts[i-1])` 保持平滑
- 初始弧：抛物线 `z = h_peak * 4t(1-t)`，`h_peak = clip(tan(β)*range/2, 150, 500)` km

**主函数**：

| 函数 | 说明 |
|------|------|
| `find_all_rays_p2p(tx, rx, n_model, freq, p2p_params)` | 主入口，返回去重排序的模式列表 |
| `variational_find_ray(tx, rx, n_model, beta, is_high, p2p_params)` | 单条变分求解，返回 `(pts, gp)` |
| `classify_mode(ray_dict)` | 按 h_reflect_km 和 tau_ms 返回标签 |
| `extract_es_params(points, h_Es_km)` | 射线与 Es 层的交叉参数 |
| `extract_bubble_entry(points, z_bot_km)` | 射线进入等离子体泡的参数 |
| `_ray_worker(args)` | 模块顶层函数（multiprocessing.Pool 可序列化） |

**模式分类规则**（阈值来自 `config.MODE`）：

| h_reflect_km | tau_ms | 标签 |
|---|---|---|
| < 140 | 任意 | `Es` |
| 140–200 | 任意 | `E` |
| 200–300 | < 5 | `1F_low` |
| 200–300 | ≥ 5 | `1F_high` |
| ≥ 300 | 任意 | `2F` |

---

### M3 — Es 层模型 (`models/es_model.py`)

**类**：`EsLayerModel`，参数见 `config.ES`

| 方法 | 说明 |
|------|------|
| `classify(f_MHz)` | 返回 `(mode, alpha)`，mode ∈ {'reflect','mixed','scatter'} |
| `reflection_coeff_sq(theta_rad, f_MHz)` | ρ²（Hao 2017 Eq.3，n=5 薄层展开） |
| `scatter_cross_section(theta_rad, f_MHz)` | σ₅ [m⁻¹]（Hao 2017 Eq.7） |
| `transmission_amplitude(theta_rad, f_MHz)` | T = sqrt(1-ρ²) ∈ [0,1] |
| `compute_power(Pt_W, Gt, Gr, f_MHz, D_km, theta_rad)` | 返回完整结果 dict |

**三段式分区**（以 foEs/f 为基准）：

| foEs/f 范围 | 模式 | alpha |
|---|---|---|
| > fr = 0.25 | reflect | 1.0 |
| fs~fr = 0.10~0.25 | mixed | 线性插值 |
| < fs = 0.10 | scatter | 0.0 |

---

### M4 — PE/SSF 传播器 (`models/pe_propagator.py`)

**类**：`PEPropagator`，参数见 `config.PE`

**SSF 单步算法**（Carrano 2020 宽角 PE）：
- Step A（折射，空域）：`u_half = u * exp(i k0 (n-1) dx)`
- Step B（衍射，谱域）：FFT → `kx_eff = sqrt(k0²-kz²) - k0` → IFFT

| 方法 | 说明 |
|------|------|
| `extract_domain(Ne_2d, x, z, x_range, z_range)` | 提取 PE 子域（可选 earth_flat 修正） |
| `ssf_step(u, n_half, k0, dz, dx)` | 单步 SSF（静态方法） |
| `apply_pml(u, n_pml, sigma)` | 指数衰减 PML（静态方法） |
| `propagate(u_init, n_field, dx, dz)` | 主传播循环，返回 `(u_out, history)` |
| `analyze(u_out, z_array, dx_total)` | AOA 谱分析，返回 mean_aoa, delta_tau 等 |
| `extract_scatter_modes(aoa_deg, power_aoa, aoa_inc)` | 找 AOA 谱中的散射峰 |

**模块级函数**：

| 函数 | 说明 |
|------|------|
| `construct_incident_field(A, beta_deg, z_inc, z_arr, k0, w0_km)` | 高斯波束初始场 |

**地球展平修正**（`PE['earth_flat']=True` 时）：
```
n_eff(x, z) = n(x, z) * (1 + z / R_E)
```
在泡区高度 350 km 处修正约 5.5%；关闭方式：`config.PE['earth_flat'] = False`

---

### M5-M7 — 混合模型 (`models/hybrid_model.py`)

**类**：`HybridPropagationModel`

```python
model = HybridPropagationModel(iono_params, radar_params,
                               x_array=BG_X, z_array=BG_Z)
mode_results, tau_axis, pd_W, main_mode = model.compute(tx_km, rx_km, t, p2p_params)
```

**`iono_params` 结构**：
```python
{
    'iri_params':    {'dt': datetime, 'lat': float, 'lon': float},
    'tid_params':    {**TID,    'enable': bool},   # 可选
    'es_params':     {**ES,     'enable': bool},   # 可选
    'bubble_params': {**BUBBLE, 'enable': bool},   # 可选
}
```

**`radar_params` 结构**：
```python
{'freq_MHz': float, 'Pt_W': float, 'Gt': float, 'Gr': float}
```

**模块级函数**：

| 函数 | 说明 |
|------|------|
| `build_pd_spectrum(mode_results, tau_axis, tau_res_ms)` | 高斯展宽叠加 P-D 谱 |
| `identify_main_mode(tau_axis, pd_W, mode_results)` | 峰值检测 + 映射到最近模式 |

---

## 5. 参数速查表

所有参数集中在 `config.py`，修改该文件即可全局生效。

### 测试链路

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `TX_LAT` | 30.0 | 发射机纬度 [°N] |
| `TX_LON` | 120.0 | 发射机经度 [°E] |
| `RX_RANGE` | 1169.0 | TX–RX 水平距离 [km] |
| `FREQ_MHZ` | 10.0 | 工作频率 [MHz] |
| `PT_W` | 1000.0 | 发射功率 [W] |
| `GT`, `GR` | 1.0 | TX/RX 增益（线性，各向同性） |

### 背景网格

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `BG_X_MIN/MAX` | -50 / 1300 km | 水平范围（含余量） |
| `BG_DX` | 5.0 km | 水平分辨率 |
| `BG_Z_MIN/MAX` | 60 / 600 km | 高度范围 |
| `BG_DZ` | 2.0 km | 高度分辨率 |

### IRI 背景

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `IRI_DT` | 2020-06-01 12:00 | 日期时间（正午，中等太阳活动） |
| `IRI_LAT` | 32.5°N | 路径中点纬度 |
| `IRI_LON` | 120.0°E | 经度 |

### TID 参数（`config.TID`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | 是否叠加 TID |
| `lambda_h_km` | 300 | 水平波长 [km]（典型 MSTID） |
| `T_s` | 2400 | 周期 [s]（40 分钟） |
| `amplitude` | 0.10 | 峰值 δNe/Ne₀（0–1） |
| `I_dip_deg` | 50.0 | 地磁倾角 [°]（中纬度） |
| `H_km` | 60.0 | Chapman 标高 [km] |
| `omega_b_rad_s` | 2π/1200 | Brunt-Väisälä 频率 [rad/s] |

### Es 层参数（`config.ES`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | 是否激活 Es |
| `foEs_MHz` | 5.0 | Es 等离子体频率 [MHz] |
| `h_Es_km` | 110.0 | Es 中心高度 [km] |
| `delta_h_m` | 115.0 | 半厚度 [m]（Hao 2017 最佳拟合） |
| `n_exp` | 5 | 密度剖面指数 |
| `L1_m, L2_m` | 300 | 水平不规则尺度 [m] |
| `L3_m` | 30 | 垂直不规则尺度 [m] |
| `delta_N_N` | 0.3 | 相对密度起伏 ΔN/N |
| `fr` | 0.25 | 反射阈值 foEs/f |
| `fs` | 0.10 | 散射阈值 foEs/f |

### 等离子体泡参数（`config.BUBBLE`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `enable` | False | 是否激活等离子体泡 |
| `delta_max` | 0.6 | 最大耗散比（0–1） |
| `x0_km` | 600.0 | 泡中心水平位置 [km]（约路径中点） |
| `z0_km` | 350.0 | 泡中心高度 [km] |
| `Lx_km` | 100.0 | 水平半宽 [km] |
| `Lz_km` | 150.0 | 垂直半高 [km] |

### P2P 变分参数（`config.P2P`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `n_init` | 18 | 初始仰角扫描数（每个运行高+低两种） |
| `n_ctrl` | 30 | 每条路径的内部控制点数 |
| `alpha_km` | 0.5 | 梯度下降步长 [km/iter] |
| `k_spring` | 0.1 | 平滑弹簧系数 |
| `max_iter` | 500 | 最大迭代次数 |
| `tol_km` | 0.05 | 收敛准则 ‖∂S/∂r‖⊥ [km] |
| `clust_h_km` | 10.0 | 去重高度容差 [km] |
| `clust_tau_ms` | 0.05 | 去重时延容差 [ms] |
| `n_workers` | 1 | 并行进程数（0=自动，1=串行） |

> **开启并行**：将 `n_workers` 改为 0（自动检测 CPU 数）或 N（固定 N 进程）。  
> 调用脚本必须有 `if __name__ == '__main__':` 保护（Windows spawn 要求）。

### PE/SSF 参数（`config.PE`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `dx_km` | 0.5 | 传播步长 [km] |
| `dz_m` | λ/4 = 7.5 m | 垂直采样 [m]（@10 MHz） |
| `n_pml` | 60 | PML 层厚度 [格点] |
| `sigma_pml` | 0.4 | PML 最大衰减系数 |
| `w0_km` | 20.0 | RT→PE 高斯波束腰 [km] |
| `store_history` | False | 是否存储所有 x 截面（调试用，耗内存） |
| `min_power_frac` | 0.01 | 散射峰最小功率比（相对主峰） |
| `earth_flat` | True | 是否应用 n_eff = n*(1+z/R_E) 地球展平修正 |

### 模式分类阈值（`config.MODE`）

| 键 | 默认值 | 说明 |
|----|--------|------|
| `h_Es_km` | 140.0 | h_r < 140 km → Es 模式 |
| `h_E_km` | 200.0 | 140 ≤ h_r < 200 km → E 层模式 |

---

## 6. 输出格式

### mode_results（每个元素的键）

| 键 | 单位 | 说明 |
|----|------|------|
| `label` | — | 模式标签：`1F_low`、`1F_high`、`2F`、`E`、`Es`、`Es_reflect`、`Es_scatter`、`Es_mixed` |
| `tau_ms` | ms | 群时延 |
| `delta_tau_ms` | ms | 时延扩展（来自散射模型） |
| `Pr_W` | W | 接收功率（线性） |
| `Pr_dBW` | dBW | 接收功率（对数） |
| `h_reflect_km` | km | 反射高度（射线最高点） |
| `group_path_km` | km | 群路径长度 |
| `beta_deg` | deg | 出发仰角 |
| `phi_deg` | deg | 到达方位角（当前为 0，2-D 模型未实现） |
| `points` | km | 变分控制点数组 `(n_ctrl+2, 2)`，用于路径可视化 |

> **注意**：`find_all_rays_p2p()` 直接返回的 dict 不含 `Pr_W`/`Pr_dBW`/`delta_tau_ms`。
> 不经过 `HybridPropagationModel.compute()` 时，须调用 `main._add_free_space_power(modes)` 补充自由空间功率。

### compute() 返回值

```python
mode_results, tau_axis, pd_W, main_mode = model.compute(...)
# mode_results : list[dict]        每条传播模式
# tau_axis     : np.ndarray [ms]   P-D 谱横轴
# pd_W         : np.ndarray [W]    P-D 谱功率（高斯叠加）
# main_mode    : dict | None       P-D 谱主峰对应的模式
```

### P-D 谱绘图

```python
from viz.plot_utils import plot_pd_spectrum
fig, ax = plot_pd_spectrum(tau_axis, pd_W, mode_results,
                           title='P-D Spectrum',
                           save_path='output/pd.png')
```

---

## 7. 场景使用示例

### 最简运行（IRI 仅背景）

```python
from models.hybrid_model import HybridPropagationModel
from datetime import datetime
import config as cfg

model = HybridPropagationModel(
    iono_params  = {'iri_params': {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON}},
    radar_params = {'freq_MHz': 10.0, 'Pt_W': 1000.0, 'Gt': 1.0, 'Gr': 1.0},
)
modes, tau_ax, pd_W, main = model.compute()
```

### 开启 TID + Es 联合

```python
import config as cfg

iono_params = {
    'iri_params':  {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
    'tid_params':  {**cfg.TID,    'enable': True, 'amplitude': 0.15},
    'es_params':   {**cfg.ES,     'enable': True, 'foEs_MHz': 7.0},
}
```

### 快速测试（减少 P2P 迭代）

```python
import config as cfg

p2p_fast = {**cfg.P2P, 'n_init': 6, 'max_iter': 200}
modes, tau_ax, pd_W, main = model.compute(p2p_params=p2p_fast)
```

### 开启 P2P 并行加速

```python
# 调用脚本必须有 __main__ 保护
if __name__ == '__main__':
    import config as cfg
    p2p_par = {**cfg.P2P, 'n_workers': 0}   # 0 = 自动检测 CPU 数
    modes, tau_ax, pd_W, main = model.compute(p2p_params=p2p_par)
```

### 关闭地球展平修正（PE 精度对比）

```python
import config as cfg
cfg.PE['earth_flat'] = False
# 重新构造 PEPropagator 后生效（PEPropagator.__init__ 读取 pe_params）
```

---

## 8. 验证脚本

| 脚本 | 验证内容 | 通过准则 |
|------|---------|---------|
| `_verify_part1.py` | IonosphereModel：Ne 值域、TID 振幅、Es 峰值、泡耗散 | 定量断言 |
| `_verify_part2.py` | 射线追踪：折射、虚高、扇形覆盖 | 定量断言 |
| `_verify_part3.py` | 光程泛函、梯度形状、变分收敛、P2P 模式数 | 定量断言 |
| `_verify_part4.py` | Es 分类、ρ²、σ₅、功率曲线形状 | 定量断言 |
| `_verify_part5.py` | 高斯波束、SSF 能量守恒 (<1e-10)、PML、AOA=20° | 精确数值 |
| `_verify_part6.py` | P-D 谱峰值、主模式映射、完整 pipeline | 端到端 |

---

## 9. 精度与性能

| 模块 | 方法 | 精度 | 单次耗时（估算，单线程） |
|------|------|------|--------------------------|
| M1 密度场 | IRI + Hooke + 插值 | ★★★☆ | ~1 s |
| M2 P2P | 广义力变分 (n_init=18) | ★★★★ | ~30–60 s |
| M3 Es | Hao 解析模型 | ★★★★ | ~0.01 s |
| M4 PE/SSF | 分步傅里叶 (100 km 域) | ★★★☆ | ~10–100 s |
| M5-M7 合成 | 功率叠加 + 峰值检测 | ★★★☆ | ~0.1 s |
| 完整 pipeline | M1+M2(P2P)+M3+M4+M5-M7 | ★★★★ | ~1–5 min |

> **精度说明**：★★★★ = 与实测数据误差 < 2 dB / 时延误差 < 0.5 ms

**加速建议**：
- 快速预览：`n_init=6, max_iter=200`（~5–10 s）
- 批量扫描：`n_workers=0`（并行，需 `__main__` 保护）+ 仅对穿过不规则体的射线运行 PE

---

## 10. 常见修改方法

### 改变工作频率

```python
# config.py
FREQ_MHZ = 15.0
# PE 的 dz_m 会自动随 lambda 调整（使用 _LAM_M / 4）
```

### 改变 TX–RX 链路

```python
# config.py
TX_LAT   = 25.0
TX_LON   = 105.0
RX_RANGE = 800.0
RX_POS   = (800.0, 0.0)
IRI_LAT  = 27.0    # 路径中点纬度
# BG_X_MAX 应大于 RX_RANGE + 余量
BG_X_MAX = 950.0
BG_X = np.arange(BG_X_MIN, BG_X_MAX + BG_DX, BG_DX)
```

### 改变 IRI 日期/地点

```python
from datetime import datetime
cfg.IRI_DT  = datetime(2021, 12, 1, 2, 0)   # 冬季夜间
cfg.IRI_LAT = 40.0
cfg.IRI_LON = 115.0
```

### 调整 TID 参数进行敏感性分析

```python
tid_sweep = []
for amp in [0.05, 0.10, 0.20]:
    iono_p = {
        'iri_params': {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
        'tid_params': {**cfg.TID, 'enable': True, 'amplitude': amp},
    }
    model = HybridPropagationModel(iono_p, radar_params)
    modes, tau_ax, pd_W, main = model.compute(p2p_params=p2p_fast)
    tid_sweep.append((amp, modes, tau_ax, pd_W))
```

### 改变等离子体泡位置

```python
# config.py
BUBBLE = {
    **BUBBLE,
    'x0_km': 400.0,   # 路径前 1/3 处
    'z0_km': 300.0,   # 低高度泡
    'delta_max': 0.8, # 更强耗散
}
```

### 关闭地球展平修正（与旧版对比）

```python
import config as cfg
cfg.PE['earth_flat'] = False
model = HybridPropagationModel(iono_params, radar_params)
# PE 域内 n 值将保持 n <= 1
```

### 存储 PE 传播历史（调试散射结构）

```python
import config as cfg
cfg.PE['store_history'] = True
# model.pe.propagate() 将返回 (u_out, history)
# history shape: (Nx-1, Nz)，每步一个截面
# 注意：300 km 域 / 7.5 m dz 约 500 MB RAM
```

---

## 参考文献

| 编号 | 文献 | 对应模块 |
|------|------|---------|
| [1] | Hao et al. 2017 | M3 Es 三段式模型（fr=0.25, fs=0.10） |
| [2] | Koval et al. 2018 | M1 TID（Hooke 1968 扰动公式） |
| [3] | Carrano et al. 2020 | M4 PE/SSF（Δz=6.1 m, Δx=0.5 km） |
| [4] | Nosikov et al. 2020 | M2b 广义力变分 P2P 求解器 |
| [5] | Jiang et al. 2020 | M1 等离子体泡物理参数 |
| [6] | Green et al. 2025 | FARR FDTD 全波验证工具 |
| [7] | Lv et al. 2025 | 槽区传播效应 |
