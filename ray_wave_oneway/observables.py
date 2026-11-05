"""Observables for exp(+i k.r - i omega t) complex transfer functions."""

import numpy as np


def phase_delay_s(plus, minus, half_step_hz):
    """Estimate delay from the principal phase difference at f +/- half_step_hz.

    The result is unambiguous only for |delay| < 1 / (4 * half_step_hz):
    for example, +/- 1 kHz permits delays below 250 microseconds.  Longer
    paths need a known reference phase removed or multistep phase unwrapping.
    """
    if half_step_hz <= 0 or abs(plus) <= 1e-15 or abs(minus) <= 1e-15:
        raise ValueError("phase delay requires nonzero responses and positive frequency step")
    return float(np.angle(plus * np.conj(minus)) / (4 * np.pi * half_step_hz))


def doppler_hz(response, slow_time_s):
    values = np.asarray(response, dtype=complex)
    times = np.asarray(slow_time_s, dtype=float)
    if (values.ndim != 1 or times.shape != values.shape or values.size < 3
            or not np.all(np.isfinite(values)) or np.any(np.abs(values) <= 1e-15)
            or not np.all(np.isfinite(times)) or not np.all(np.diff(times) > 0)):
        raise ValueError("Doppler requires nonzero complex responses at increasing times")
    return -np.gradient(np.unwrap(np.angle(values)), times, edge_order=2) / (2 * np.pi)


def power_delay(frequency_hz, response):
    frequency = np.asarray(frequency_hz, dtype=float)
    values = np.asarray(response, dtype=complex)
    if (frequency.ndim != 1 or values.shape != frequency.shape or frequency.size < 2
            or not np.all(np.isfinite(frequency)) or not np.all(np.isfinite(values))):
        raise ValueError("power-delay input must contain finite frequency and response vectors")
    spacing = np.diff(frequency)
    if spacing[0] <= 0 or not np.allclose(spacing, spacing[0], rtol=1e-10, atol=1e-8):
        raise ValueError("power-delay frequencies must be uniformly increasing")
    delay = np.fft.fftshift(np.fft.fftfreq(values.size, d=spacing[0]))
    pulse = np.fft.fftshift(np.fft.fft(values, norm="ortho")) / np.sqrt(values.size)
    return delay, np.abs(pulse)**2


def power_doppler(slow_time_s, response):
    times = np.asarray(slow_time_s, dtype=float)
    values = np.asarray(response, dtype=complex)
    if (times.ndim != 1 or values.shape != times.shape or times.size < 2
            or not np.all(np.isfinite(times)) or not np.all(np.isfinite(values))):
        raise ValueError("power-Doppler input must contain finite time and response vectors")
    spacing = np.diff(times)
    if spacing[0] <= 0 or not np.allclose(spacing, spacing[0], rtol=1e-10, atol=1e-8):
        raise ValueError("slow times must be uniformly increasing")
    frequency = np.fft.fftshift(np.fft.fftfreq(values.size, d=spacing[0]))
    spectrum = np.fft.fftshift(np.fft.ifft(values, norm="ortho")) / np.sqrt(values.size)
    return frequency, np.abs(spectrum)**2
