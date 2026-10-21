"""Isotope and diatomic reduced masses.

Isotope data are provided by :mod:`periodictable` in unified atomic mass
units.  Nuclear masses are approximated from neutral-atom masses by removing
the rest mass of every electron.  Electronic binding energies are not
included in that conversion.
"""

from __future__ import annotations

import re
from collections.abc import Sequence

import periodictable

from .scipy_codata import get_codata_table


_ISOTOPE_PATTERNS = (
    re.compile(r"^(?P<symbol>[A-Za-z]{1,3})-(?P<mass_number>\d+)$"),
    re.compile(r"^(?P<mass_number>\d+)-(?P<symbol>[A-Za-z]{1,3})$"),
    re.compile(r"^(?P<mass_number>\d+)(?P<symbol>[A-Za-z]{1,3})$"),
)
_HYDROGEN_ALIASES = {"D": ("H", 2), "T": ("H", 3)}
_UNIT_ALIASES = {
    "u": "u",
    "da": "u",
    "dalton": "u",
    "daltons": "u",
    "electron_mass": "electron_mass",
    "electron_masses": "electron_mass",
    "au": "electron_mass",
    "a.u.": "electron_mass",
}


def _parse_isotope(isotope: str):
    if not isinstance(isotope, str):
        raise TypeError("isotope must be a string such as 'Li-7' or 'H-1'.")

    specification = isotope.strip()
    if specification in _HYDROGEN_ALIASES:
        symbol, mass_number = _HYDROGEN_ALIASES[specification]
    else:
        match = next(
            (
                match
                for pattern in _ISOTOPE_PATTERNS
                if (match := pattern.fullmatch(specification)) is not None
            ),
            None,
        )
        if match is None:
            raise ValueError(
                f"Invalid isotope specification '{isotope}'. Use an explicit "
                "isotope such as 'Li-7' or 'H-1'."
            )
        symbol = match.group("symbol").capitalize()
        mass_number = int(match.group("mass_number"))

    try:
        element = getattr(periodictable, symbol)
        result = element[mass_number]
    except (AttributeError, KeyError, ValueError) as exc:
        raise ValueError(f"Unknown isotope '{isotope}'.") from exc

    if result.mass is None:
        raise ValueError(f"No mass is available for isotope '{isotope}'.")
    return result


def _mass_constants(codata_set: str) -> tuple[float, float]:
    constants, _ = get_codata_table(codata_set)
    electron_mass_kg = constants["electron mass"][0]
    atomic_mass_constant_kg = constants["atomic mass constant"][0]
    electron_mass_u = electron_mass_kg / atomic_mass_constant_kg
    return electron_mass_u, atomic_mass_constant_kg / electron_mass_kg


def _normalize_units(units: str) -> str:
    if not isinstance(units, str):
        raise TypeError("units must be a string.")
    try:
        return _UNIT_ALIASES[units.lower()]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported mass units '{units}'. Use 'u' or 'electron_mass'."
        ) from exc


def isotope_mass(
    isotope: str,
    *,
    units: str = "electron_mass",
    mass_kind: str = "atomic",
    codata_set: str = "current",
) -> float:
    """Return an isotope mass.

    Parameters
    ----------
    isotope
        Explicit isotope specification, for example ``"Li-7"``, ``"7-Li"``,
        ``"H-1"``, or the hydrogen aliases ``"D"`` and ``"T"``.
    units
        ``"u"`` for unified atomic mass units or ``"electron_mass"`` for
        atomic units of mass. Common aliases such as ``"Da"`` and ``"au"``
        are accepted.
    mass_kind
        ``"atomic"`` returns the tabulated neutral-atom mass. ``"nuclear"``
        subtracts the rest masses of all electrons. This nuclear conversion
        neglects electronic binding energy.
    codata_set
        CODATA set used to convert between unified atomic mass units and
        electron masses. See :func:`matter_field_params.get_codata_table`.
    """
    isotope_data = _parse_isotope(isotope)
    electron_mass_u, u_to_electron_mass = _mass_constants(codata_set)

    normalized_kind = mass_kind.lower()
    if normalized_kind == "atomic":
        mass_u = isotope_data.mass
    elif normalized_kind == "nuclear":
        mass_u = (
            isotope_data.mass
            - isotope_data.element.number * electron_mass_u
        )
    else:
        raise ValueError(
            f"Unsupported mass_kind '{mass_kind}'. Use 'atomic' or 'nuclear'."
        )

    normalized_units = _normalize_units(units)
    if normalized_units == "u":
        return float(mass_u)
    return float(mass_u * u_to_electron_mass)


def nuclear_mass(
    isotope: str,
    *,
    units: str = "electron_mass",
    codata_set: str = "current",
) -> float:
    """Return an approximate bare-nucleus mass for an isotope.

    The approximation subtracts all electron rest masses from the tabulated
    neutral-atom mass and neglects electronic binding energy.
    """
    return isotope_mass(
        isotope,
        units=units,
        mass_kind="nuclear",
        codata_set=codata_set,
    )


def reduced_mass(
    isotopes: Sequence[str],
    *,
    units: str = "electron_mass",
    mass_kind: str = "nuclear",
    codata_set: str = "current",
) -> float:
    """Return the reduced mass of two explicitly specified isotopes."""
    if isinstance(isotopes, str) or len(isotopes) != 2:
        raise ValueError("isotopes must contain exactly two isotope strings.")

    mass_1, mass_2 = (
        isotope_mass(
            isotope,
            units=units,
            mass_kind=mass_kind,
            codata_set=codata_set,
        )
        for isotope in isotopes
    )
    return mass_1 * mass_2 / (mass_1 + mass_2)
