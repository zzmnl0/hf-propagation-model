"""Parameters for the isolated controlled one-way scattering case."""

from dataclasses import dataclass
from pathlib import Path


# These values belong to this route; they are not read from the original AH
# UPDATE2 profile.  Pass the finite-difference values explicitly to each model.
ONEWAY_L2_SETTINGS = {
    "rtol": 1.0e-10,
    "atol_r_km": 1.0e-8,
    "atol_p": 1.0e-12,
    "ds_min_km": 0.002,
    "ds_max_km": 0.5,
}
ONEWAY_AH_FD_SETTINGS = {
    "spatial_fd_base_km": 0.03125,
    "spatial_fd_scales": (0.5, 1.0, 2.0),
    "spatial_fd_stability_rtol": 0.02,
    "spatial_fd_stability_atol_n2_per_km": 2.0e-5,
}


@dataclass(frozen=True)
class ControlledCase:
    frequency_hz: float = 10e6
    plasma_coefficient: float = 80.61638587963631
    earth_radius_km: float = 6371.0
    background_X_at_center: float = .2
    magnetic_field_T: float = 45e-6
    medium_time_s: float = 37.0
    slab_height_km: float = 250.0
    period_u_m: float = 500.0
    period_v_m: float = 600.0
    before_m: float = 100.0
    after_m: float = 100.0
    irregularity_fraction: float = .05
    grid_size: int = 5
    sweep_half_span_hz: float = 500e3
    sweep_count: int = 9
    slow_time_step_s: float = .25
    slow_time_count: int = 9

    @property
    def background_Ne_m3(self):
        return self.background_X_at_center * self.frequency_hz**2 / self.plasma_coefficient

    @property
    def output_dir(self):
        return Path(__file__).resolve().parents[1] / "output2" / "ray_wave_oneway" / "s4" / "observables"
