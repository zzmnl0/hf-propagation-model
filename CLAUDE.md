# HF Shortwave Propagation Hybrid Model

Research project 李慧敏-22: 复杂电离层环境短波反射/散射混合传播机理及建模研究

## Directory structure

```
code/
├── config.py                  # All physical parameters and grid settings
├── utils.py                   # Physics utility functions (EM, power, coords)
├── main.py                    # Entry point; uncomment scenarios to run
├── models/                    # Physics modules M1-M7
│   ├── ionosphere_model.py    # M1: IRI + TID / Es / plasma bubble  [DONE]
│   ├── ray_tracer.py          # M2: Haselgrove RK4 ray tracer        [Part 2]
│   ├── point_to_point.py      # M2: P2P variational solver (Nosikov)  [Part 3]
│   ├── es_model.py            # M3: Es three-segment model (Hao 2017) [Part 4]
│   ├── pe_propagator.py       # M4: PE/SSF bubble scatter (Carrano)   [Part 5]
│   └── hybrid_model.py        # M5-M7: full pipeline                  [Part 6]
├── viz/                       # Visualisation
│   ├── plot_utils.py          # Shared plot helpers (ray fan, PD, Ne field)
│   └── plot_ne_background.py  # 4-panel Ne background figures -> output/
├── tests/                     # Verification scripts
│   ├── _verify_part0.py       # Import / config sanity check
│   ├── _verify_part1.py       # IonosphereModel full verification
│   └── _probe_iri2016.py      # iri2016 API probe
└── output/                    # Generated PNG files (gitignored)
```

## Running

```
# Environment verification
conda run -n pytorch_cpu python main.py

# Part 1 verification (IRI + TID + Es + bubble)
conda run -n pytorch_cpu python tests/_verify_part1.py

# Background Ne plots -> output/ne_*.png
conda run -n pytorch_cpu python viz/plot_ne_background.py
```

## Key conventions

- All distances in **km**, heights in **km**, frequency in **MHz**, time in **s**
- 2-D grid: `x` = horizontal distance from TX [km], `z` = height [km]
- Ne in **m^-3**; refractive index `n` is real (isotropic, no magnetic field)
- IRI via `iri2016.IRI(dt, (z_min, z_max, dz), lat, lon)['ne'].values`  (NOT iricore)
- Inter-module imports inside `models/` use relative form: `from .ray_tracer import`
- Visualization functions live in `viz/plot_utils.py`, not `utils.py`
- Output PNGs go to `output/`

## Implementation status

| Part | Module | Status |
|------|--------|--------|
| 1 | `models/ionosphere_model.py` | **Done** |
| 2 | `models/ray_tracer.py` | **Done** |
| 3 | `models/point_to_point.py` | **Done** |
| 4 | `models/es_model.py` | **Done** |
| 5 | `models/pe_propagator.py` | Stub |
| 6 | `models/hybrid_model.py` | Stub |
