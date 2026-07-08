# stage1 P0 traceability matrix

| ID | Requirement | Code | Verification | Status |
| --- | --- | --- | --- | --- |
| P0-1 | Default feature gates remain off | config.py | tests/_verify_update2_all.py | implemented |
| P0-2 | 3-D NeProvider and parameterized ionosphere | models/ne_provider.py; models/parameterized_ionosphere_3d.py | tests/_verify_ne_provider_3d.py | implemented |
| P0-3 | AH3D Hamiltonian tracing diagnostics | models/ray_tracer_ah_3d.py | tests/_verify_AH3D_physics.py | implemented |
| P0-4 | Reproducible compact 3-D parameter sweep | analysis/parameter_sweep_3d.py | tests/_verify_sweep_3d_reproducibility.py | implemented |
