"""P2-2 verification: single-station OTH radar clutter physics.

Checks:
  1. sigma0=0 gives zero echo power.
  2. zero radial velocity gives zero Doppler.
  3. two-way delay equals 2x one-way delay.
  4. received power follows R^-4 monostatic radar scaling.
  5. default RADAR_CLUTTER feature gate is off.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from legacy.models.radar_clutter import RadarClutterModel, apply_clutter_to_mode


def _mode(group_path_km=1000.0, tau_ms=3.0, area_km2=50.0):
    return {
        "label": "1F",
        "tau_ms": tau_ms,
        "group_path_km": group_path_km,
        "A_tube_km2": area_km2,
        "beta_deg": 30.0,
        "Pr_W": 1.0,
        "Pr_dBW": 0.0,
    }


def run():
    n_fail = 0
    print("=" * 55)
    print("  P2-2 verification: radar clutter physics")
    print("=" * 55)

    model = RadarClutterModel(freq_MHz=cfg.FREQ_MHZ)

    print("\nCheck 1: sigma0=0 -> Pr=0")
    r0 = model.evaluate_mode(_mode(), params={"sigma0_dB": -400.0})
    ok1 = r0["clutter_Pr_W"] < 1e-35 and r0["sigma_eff_m2"] < 1e-20
    print("  Pr={:.3e} W sigma_eff={:.3e} m2 {}".format(
        r0["clutter_Pr_W"], r0["sigma_eff_m2"], "OK" if ok1 else "FAIL"))
    if not ok1:
        n_fail += 1

    print("\nCheck 2: v_radial=0 -> Doppler=0")
    r1 = model.evaluate_mode(_mode(), params={"radial_velocity_ms": 0.0})
    ok2 = abs(r1["doppler_hz"]) < 1e-12
    print("  Doppler={:.6f} Hz {}".format(r1["doppler_hz"], "OK" if ok2 else "FAIL"))
    if not ok2:
        n_fail += 1

    print("\nCheck 3: tau_2way = 2*tau_oneway")
    r2 = model.evaluate_mode(_mode(tau_ms=4.25))
    ok3 = abs(r2["tau_2way_ms"] - 8.5) < 1e-12
    print("  tau1={:.3f} ms tau2={:.3f} ms {}".format(
        r2["tau_oneway_ms"], r2["tau_2way_ms"], "OK" if ok3 else "FAIL"))
    if not ok3:
        n_fail += 1

    print("\nCheck 4: monostatic power follows R^-4")
    near = model.evaluate_mode(_mode(group_path_km=1000.0), params={"sigma0_dB": -20.0})
    far = model.evaluate_mode(_mode(group_path_km=2000.0), params={"sigma0_dB": -20.0})
    ratio = near["clutter_Pr_W"] / max(far["clutter_Pr_W"], 1e-300)
    ok4 = abs(ratio - 16.0) / 16.0 < 1e-12
    print("  Pr(1000km)/Pr(2000km)={:.6f} {}".format(ratio, "OK" if ok4 else "FAIL"))
    if not ok4:
        n_fail += 1

    print("\nCheck 5: feature gate default off and apply_to_power explicit")
    gate = cfg.UPDATE2.get("RADAR_CLUTTER", {})
    d0 = apply_clutter_to_mode(
        _mode(),
        freq_MHz=cfg.FREQ_MHZ,
        params={"sigma0_dB": -20.0},
        apply_to_power=False,
    )
    d1 = apply_clutter_to_mode(
        _mode(),
        freq_MHz=cfg.FREQ_MHZ,
        params={"sigma0_dB": -20.0},
        apply_to_power=True,
    )
    ok5 = (
        gate.get("enable") is False
        and gate.get("apply_to_power") is False
        and abs(d0["Pr_W"] - 1.0) < 1e-12
        and d1["Pr_W"] == d1["clutter_Pr_W"]
    )
    print("  enable={} apply_to_power={} {}".format(
        gate.get("enable"), gate.get("apply_to_power"), "OK" if ok5 else "FAIL"))
    if not ok5:
        n_fail += 1

    print("\n{} fail(s)\n".format(n_fail))
    return n_fail


if __name__ == "__main__":
    sys.exit(run())
