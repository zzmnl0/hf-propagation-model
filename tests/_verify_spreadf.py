"""
Phase 3 verification: spread-F power-law phase screen (Rino 1979).

Checks:
  1. Cs=0 -> Ne_2d unchanged after SpreadFModel.apply()
  2. Cs>0 -> Ne perturbed at h_screen_km, near-zero far from screen
  3. FFT power spectrum slope ~ -(p+1) in mid-k range
  4. IonosphereModel with spread_f enabled produces different Ne than without
  5. Larger Cs -> larger Ne perturbation (monotone)

Run:
    conda run -n pytorch_cpu python tests/_verify_spreadf.py
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from models.spread_f_model import SpreadFModel
from models.ionosphere_model import IonosphereModel
import config as cfg


def run():
    n_fail = 0
    x_km   = cfg.BG_X
    z_km   = cfg.BG_Z

    # Build background Ne once
    iono   = IonosphereModel()
    Ne_bg, _ = iono.build_Ne_field(x_km, z_km)

    # ── Check 1: Cs=0 -> no perturbation ─────────────────────────────────────
    print("[1] Cs=0 -> Ne_2d unchanged")
    sfm0  = SpreadFModel(Cs=0.0, seed=42)
    Ne_0  = sfm0.apply(Ne_bg.copy(), x_km, z_km)
    diff0 = float(np.max(np.abs(Ne_0 - Ne_bg)))
    ok    = diff0 == 0.0
    print("    max |dNe| = {:.2e}  {}".format(diff0, 'OK' if ok else 'FAIL'))
    if not ok:
        n_fail += 1

    # ── Check 2: Cs>0 -> Ne perturbed near h_screen, small far from it ───────
    print("[2] Cs>0 -> perturbation near h_screen_km=300 km")
    sfm1   = SpreadFModel(Cs=1e-3, p=3.0, h_screen_km=300.0, L0_km=50.0, seed=42)
    Ne_1   = sfm1.apply(Ne_bg.copy(), x_km, z_km)

    iz300  = int(np.argmin(np.abs(z_km - 300.0)))
    iz_far = int(np.argmin(np.abs(z_km - 150.0)))   # far from screen (150 km)
    dNe_at = float(np.mean(np.abs(Ne_1[:, iz300] - Ne_bg[:, iz300])))
    dNe_far= float(np.mean(np.abs(Ne_1[:, iz_far] - Ne_bg[:, iz_far])))
    ok = dNe_at > 0 and dNe_at > dNe_far * 5
    print("    mean|dNe| at 300km = {:.2e}   at 150km = {:.2e}   {}".format(
        dNe_at, dNe_far, 'OK' if ok else 'FAIL'))
    if not ok:
        n_fail += 1

    # ── Check 3: FFT power spectrum slope ~ -(p+1) ───────────────────────────
    print("[3] Power spectrum slope ~ -(p+1) = -4.0 for p=3")
    phi = sfm1._gen_screen(x_km)
    N   = len(phi)
    dx  = float(x_km[1] - x_km[0])
    k_cyc = np.fft.rfftfreq(N, d=dx)          # cycles/km
    P     = np.abs(np.fft.rfft(phi)) ** 2

    k_lo = 1.0 / sfm1.L0_km                   # outer-scale frequency (cycles/km)
    k_hi = 0.05                                # upper limit to avoid Nyquist noise
    mask = (k_cyc > k_lo) & (k_cyc < k_hi) & (P > 0)
    if mask.sum() >= 5:
        slope = float(np.polyfit(np.log10(k_cyc[mask]),
                                 np.log10(P[mask]), 1)[0])
        expected = -(sfm1.p + 1)
        ok = abs(slope - expected) < 1.0
        print("    slope = {:.2f}   expected ~ {:.1f}   diff = {:.2f}   {}".format(
            slope, expected, abs(slope - expected), 'OK' if ok else 'FAIL'))
        if not ok:
            n_fail += 1
    else:
        print("    WARN: too few points in k-range for slope fit (mask.sum={})".format(
            mask.sum()))

    # ── Check 4: IonosphereModel with spread_f -> different Ne ───────────────
    print("[4] IonosphereModel with spread_f enabled produces different Ne")
    iono_sf  = IonosphereModel(
        spread_f_params={**cfg.SPREAD_F, 'enable': True, 'seed': 42})
    Ne_sf, _ = iono_sf.build_Ne_field(x_km, z_km)
    max_diff = float(np.max(np.abs(Ne_sf - Ne_bg)))
    ok = max_diff > 0
    print("    max |dNe| = {:.2e}  {}".format(max_diff, 'OK' if ok else 'FAIL'))
    if not ok:
        n_fail += 1

    # ── Check 5: larger Cs -> larger perturbation (monotone) ─────────────────
    print("[5] Larger Cs -> larger Ne perturbation")
    dNe_vals = []
    for Cs in (1e-4, 5e-4, 1e-3, 5e-3):
        sfm_c = SpreadFModel(Cs=Cs, p=3.0, h_screen_km=300.0,
                             L0_km=50.0, seed=42)
        Ne_c = sfm_c.apply(Ne_bg.copy(), x_km, z_km)
        dNe_vals.append(float(np.mean(np.abs(Ne_c[:, iz300] - Ne_bg[:, iz300]))))
        print("    Cs={:.0e}  mean|dNe| at 300km = {:.2e}".format(Cs, dNe_vals[-1]))
    ok = all(dNe_vals[i] < dNe_vals[i+1] for i in range(len(dNe_vals)-1))
    print("    Monotone increasing  {}".format('OK' if ok else 'FAIL'))
    if not ok:
        n_fail += 1

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n" + "=" * 50)
    if n_fail == 0:
        print("Phase 3 spread-F verification PASSED.")
    else:
        print("{} check(s) FAILED.".format(n_fail))
    return n_fail


if __name__ == '__main__':
    sys.exit(run())
