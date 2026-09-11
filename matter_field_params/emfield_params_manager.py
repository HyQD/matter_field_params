"""Helpers for electromagnetic-field parameter conversions.

In the dipole approximation, a pulse is commonly represented as
E(t) = E0 * f(t) * sin(omega * t + phi).

This module provides minimal conversions between common lab units and atomic
units for the carrier parameters E0 and omega.
"""

from __future__ import annotations

import numpy as np
from .scipy_codata import get_codata_table


def _constant_value(table, *keys: str) -> float:
    for key in keys:
        if key in table:
            return table[key][0]
    tried = ", ".join(keys)
    raise KeyError(f"None of these constants were found in CODATA table: {tried}")


def _emfield_units_from_codata(codata_set: str):
    table, resolved_codata_set = get_codata_table(codata_set)

    atomic_unit_of_time_s = _constant_value(table, "atomic unit of time")
    atomic_unit_of_efield_v_per_m = _constant_value(table, "atomic unit of electric field")
    speed_of_light_m_per_s = _constant_value(table, "speed of light in vacuum")
    vacuum_permittivity = _constant_value(
        table,
        "vacuum electric permittivity",
        "electric constant",
    )

    atomic_unit_of_intensity_w_per_m2 = (
        0.5
        * vacuum_permittivity
        * speed_of_light_m_per_s
        * atomic_unit_of_efield_v_per_m**2
    )
    atomic_unit_of_intensity_w_per_cm2 = atomic_unit_of_intensity_w_per_m2 / 1.0e4

    return {
        "codata_set": resolved_codata_set,
        "atomic_unit_of_time_s": atomic_unit_of_time_s,
        "atomic_unit_of_intensity_w_per_cm2": atomic_unit_of_intensity_w_per_cm2,
        "speed_of_light_m_per_s": speed_of_light_m_per_s,
    }


def _to_numpy_float(x):
    return np.asarray(x, dtype=float)


def _return_like_input(original, value):
    if np.isscalar(original):
        return float(np.asarray(value))
    return value


def intensity_wcm2_to_au(intensity_wcm2, codata_set: str = "current"):
    """Convert intensity from W/cm^2 to atomic units.

    Parameters
    ----------
    intensity_wcm2
        Scalar or array intensity in W/cm^2.
    codata_set : str
        CODATA set understood by get_codata_table.
    """
    intensity_wcm2_arr = _to_numpy_float(intensity_wcm2)
    if np.any(intensity_wcm2_arr < 0.0):
        raise ValueError("Intensity must be non-negative.")

    units = _emfield_units_from_codata(codata_set)
    intensity_au = intensity_wcm2_arr / units["atomic_unit_of_intensity_w_per_cm2"]
    return _return_like_input(intensity_wcm2, intensity_au)


def intensity_wcm2_to_field_au(intensity_wcm2, codata_set: str = "current"):
    """Convert intensity in W/cm^2 to electric-field amplitude E0 in a.u."""
    intensity_au = _to_numpy_float(
        intensity_wcm2_to_au(intensity_wcm2, codata_set=codata_set)
    )
    field_au = np.sqrt(2.0 * intensity_au)
    return _return_like_input(intensity_wcm2, field_au)


def field_au_to_intensity_wcm2(field_au, codata_set: str = "current"):
    """Convert electric-field amplitude E0 in a.u. to intensity in W/cm^2."""
    field_au_arr = _to_numpy_float(field_au)
    intensity_au = 0.5 * np.square(field_au_arr)
    units = _emfield_units_from_codata(codata_set)
    intensity_wcm2 = intensity_au * units["atomic_unit_of_intensity_w_per_cm2"]
    return _return_like_input(field_au, intensity_wcm2)


def wavelength_nm_to_omega_au(wavelength_nm, codata_set: str = "current"):
    """Convert carrier wavelength in nm to angular frequency omega in a.u."""
    wavelength_nm_arr = _to_numpy_float(wavelength_nm)
    if np.any(wavelength_nm_arr <= 0.0):
        raise ValueError("Wavelength must be strictly positive.")

    units = _emfield_units_from_codata(codata_set)
    wavelength_m = wavelength_nm_arr * 1.0e-9
    omega_si = 2.0 * np.pi * units["speed_of_light_m_per_s"] / wavelength_m
    omega_au = omega_si * units["atomic_unit_of_time_s"]
    return _return_like_input(wavelength_nm, omega_au)
