"""Appleton-Hartree dispersion roots and Hamiltonian derivatives.

The real dispersion surface determines the ray.  Collisional imaginary parts
belong to the separate attenuation integral; this module does not create
complex-coordinate rays.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from config import K_FP

C_KM_S = 299792.458


@dataclass(frozen=True)
class AHRoot:
    n2: complex
    branch_id: str
    status: str
    branch_gap: float
    branch_continuity: float


def _unit_analytic(v) -> np.ndarray:
    a = np.asarray(v)
    norm = np.sqrt(np.dot(a, a))
    if abs(norm) < 1e-30:
        raise ValueError("k vector must be non-zero")
    return a / norm


def appleton_hartree_roots(
    omega, Ne_m3: float, k_vec, b_hat, gyro_MHz: float, collision_s: float = 0.0,
    *, plasma_coefficient: float | None = None,
) -> dict[str, complex]:
    """Return the two cold-plasma magnetoionic phase-index-squared roots."""
    f_hz = omega / (2.0 * np.pi)
    coefficient = K_FP**2 if plasma_coefficient is None else float(plasma_coefficient)
    if not np.isfinite(coefficient) or coefficient <= 0.0:
        raise ValueError("plasma_coefficient must be finite and positive")
    X = coefficient * max(float(Ne_m3), 0.0) / (f_hz * f_hz)
    Y = float(gyro_MHz) * 1e6 / f_hz
    Z = float(collision_s) / omega
    U = 1.0 - 1j * Z
    kh = _unit_analytic(k_vec)
    bh = np.asarray(b_hat, dtype=float)
    cos_a = np.dot(kh, bh)
    yt2 = Y * Y * (1.0 - cos_a * cos_a)
    yl2 = Y * Y * cos_a * cos_a
    ux = U - X
    if abs(ux) < 1e-14:
        return {"O": complex(np.nan), "X": complex(np.nan)}
    half = yt2 / (2.0 * ux)
    disc = np.sqrt(yt2 * yt2 / (4.0 * ux * ux) + yl2 + 0j)
    return {
        "O": 1.0 - X / (U - half + disc),
        "X": 1.0 - X / (U - half - disc),
    }


def select_root(roots: dict[str, complex], branch: str, previous_n2=None,
                near_cutoff_n2: float = 2e-2, resonance_limit: float = 1e6) -> AHRoot:
    """Select a continuous root and classify its physical propagation state."""
    requested = str(branch).upper()
    if requested not in roots:
        raise ValueError("branch must be 'O' or 'X'")
    finite = {key: value for key, value in roots.items()
              if np.isfinite(value.real) and np.isfinite(value.imag)}
    if not finite:
        return AHRoot(complex(np.nan), requested, "invalid_dispersion", np.nan, np.inf)
    selected_id = requested
    if previous_n2 is not None:
        selected_id = min(
            finite,
            key=lambda key: (abs(finite[key] - previous_n2), key != requested),
        )
    n2 = finite.get(selected_id, complex(np.nan))
    other = "X" if selected_id == "O" else "O"
    gap = abs(n2 - roots.get(other, complex(np.nan)))
    continuity = 0.0 if previous_n2 is None else abs(n2 - previous_n2)
    if selected_id != requested:
        status = "branch_lost"
    elif abs(n2) > resonance_limit:
        status = "resonance_excluded"
    elif n2.real <= 0.0:
        status = "evanescent"
    elif n2.real <= near_cutoff_n2:
        status = "near_cutoff"
    else:
        status = "ok"
    return AHRoot(n2, selected_id, status, float(gap), float(continuity))


def dispersion_D(k_vec, omega, n2: complex) -> complex:
    """Haselgrove dispersion function D=|k|^2-(omega/c)^2 n^2."""
    k = np.asarray(k_vec)
    return np.dot(k, k) - (omega / C_KM_S) ** 2 * n2


def dispersion_partials(root_fn, k_vec, omega: float) -> tuple[np.ndarray, float]:
    """Compute D_k and D_omega by complex-step differentiation of one root."""
    k = np.asarray(k_vec, dtype=float)
    hk = max(1e-7, 1e-7 * float(np.linalg.norm(k)))
    hw = max(1e-2, 1e-8 * abs(float(omega)))
    # Symmetric complex steps avoid contamination by a physical imaginary part.
    dk = np.empty(3, dtype=float)
    for i in range(3):
        kp = k.astype(complex)
        km = k.astype(complex)
        kp[i] += 1j * hk
        km[i] -= 1j * hk
        dk[i] = ((root_fn(kp, omega) - root_fn(km, omega)) / (2j * hk)).real
    dw = ((root_fn(k, omega + 1j * hw) - root_fn(k, omega - 1j * hw)) / (2j * hw)).real
    return dk, float(dw)


def group_velocity(root_fn, k_vec, omega: float) -> np.ndarray:
    dk, dw = dispersion_partials(root_fn, k_vec, omega)
    if not np.isfinite(dw) or abs(dw) < 1e-20:
        return np.full(3, np.nan)
    return -dk / dw
