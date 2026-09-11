"""Top-level package for HyQD-matter-field-params."""

from .emfield_params_manager import (
    field_au_to_intensity_wcm2,
    intensity_wcm2_to_au,
    intensity_wcm2_to_field_au,
    wavelength_nm_to_omega_au,
)
from .particle_params_manager import (
    get_reduced_masses_and_charges_3body,
    predefined_3body_systems,
)
from .scipy_codata import get_codata_table

__all__ = [
    "field_au_to_intensity_wcm2",
    "get_codata_table",
    "get_reduced_masses_and_charges_3body",
    "intensity_wcm2_to_au",
    "intensity_wcm2_to_field_au",
    "predefined_3body_systems",
    "wavelength_nm_to_omega_au",
]

__version__ = "0.0.0"
