import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import config as cfg
from utils import free_space_loss_dB, to_dBW
from models.hybrid_model import HybridPropagationModel
import models.hybrid_model as hybrid_model


def ok(label, cond):
    print("  {} {}".format(label, "OK" if cond else "FAIL"))
    return 0 if cond else 1


class NoPE:
    def extract_domain(self, *args, **kwargs):
        raise AssertionError("PE extract_domain called while disabled")
    def propagate(self, *args, **kwargs):
        raise AssertionError("PE propagate called while disabled")
    def analyze(self, *args, **kwargs):
        raise AssertionError("PE analyze called while disabled")


class FlatIono:
    def __init__(self, Ne):
        self.Ne = Ne
    def build_Ne_field(self, x_array, z_array, t=0.0):
        return self.Ne.copy(), {}


def run():
    print("=" * 60)
    print("  P0 verification: PE_COUPLING disabled semantics")
    print("=" * 60)
    n_fail = 0
    old_find = hybrid_model.find_all_rays_p2p
    old_cfg = dict(cfg.UPDATE2["PE_COUPLING"])
    try:
        cfg.UPDATE2["PE_COUPLING"].update({"enable": False, "apply_to_power": False, "mode": "local_2p5d"})
        x = np.linspace(0.0, 100.0, 11)
        z = np.linspace(0.0, 300.0, 31)
        Ne = np.zeros((len(x), len(z)), dtype=float)
        pts = np.array([[0.0, 0.0], [20.0, 220.0], [50.0, 220.0], [100.0, 0.0]])
        gp_km = 100.0

        def fake_find(*args, **kwargs):
            return [{"group_path_km": gp_km, "tau_ms": gp_km / cfg.C_KMS, "h_reflect_km": 220.0, "label": "fake_F", "points": pts, "beta_deg": 30.0}]

        hybrid_model.find_all_rays_p2p = fake_find
        model = HybridPropagationModel(
            {"bubble_params": {"enable": True, "z0_km": 230.0, "Lz_km": 30.0, "Lx_km": 20.0}},
            {"Pt_W": 1.0, "Gt": 1.0, "Gr": 1.0, "freq_MHz": cfg.FREQ_MHZ},
            x_array=x,
            z_array=z,
        )
        model.iono = FlatIono(Ne)
        model.pe = NoPE()
        modes, _, _, _ = model.compute(tx_km=(0.0, 0.0), rx_km=(100.0, 0.0))
        expected = 10.0 ** (-free_space_loss_dB(gp_km, cfg.FREQ_MHZ) / 10.0)
        n_fail += ok("one fake mode returned:", len(modes) == 1)
        n_fail += ok("Pr_W unchanged by PE:", abs(modes[0]["Pr_W"] - expected) <= max(1e-30, expected * 1e-12))
        n_fail += ok("Pr_dBW consistent:", abs(modes[0]["Pr_dBW"] - to_dBW(expected)) < 1e-9)
        n_fail += ok("no pe_* diagnostics when disabled:", not any(k.startswith("pe_") for k in modes[0]))
        n_fail += ok("delta_tau_ms unchanged:", abs(float(modes[0]["delta_tau_ms"])) < 1e-15)
    finally:
        hybrid_model.find_all_rays_p2p = old_find
        cfg.UPDATE2["PE_COUPLING"].clear()
        cfg.UPDATE2["PE_COUPLING"].update(old_cfg)
    print("\n{} fail(s)\n".format(n_fail))
    return n_fail


if __name__ == "__main__":
    sys.exit(run())
