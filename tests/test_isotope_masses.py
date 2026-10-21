import periodictable
import pytest

from matter_field_params import isotope_mass, nuclear_mass, reduced_mass
from matter_field_params.scipy_codata import get_codata_table


def test_atomic_isotope_mass_in_u():
    assert isotope_mass("Li-7", units="u") == pytest.approx(
        periodictable.Li[7].mass
    )


@pytest.mark.parametrize("isotope", ["Li-7", "7-Li", "7Li"])
def test_supported_isotope_notation(isotope):
    assert isotope_mass(isotope, units="u") == pytest.approx(
        periodictable.Li[7].mass
    )


def test_hydrogen_isotope_alias():
    assert isotope_mass("D", units="u") == pytest.approx(
        periodictable.H[2].mass
    )


def test_nuclear_mass_subtracts_all_electron_rest_masses():
    constants, _ = get_codata_table("current")
    electron_mass_u = (
        constants["electron mass"][0]
        / constants["atomic mass constant"][0]
    )
    expected = periodictable.Li[7].mass - 3 * electron_mass_u

    assert nuclear_mass("Li-7", units="u") == pytest.approx(expected)


def test_reduced_mass_uses_selected_mass_kind_and_units():
    mass_li = isotope_mass("Li-7", mass_kind="nuclear")
    mass_h = isotope_mass("H-1", mass_kind="nuclear")

    assert reduced_mass(("Li-7", "H-1")) == pytest.approx(
        mass_li * mass_h / (mass_li + mass_h)
    )


@pytest.mark.parametrize("isotope", ["Li", "Li-seven", "unknown-7"])
def test_invalid_isotope_specification(isotope):
    with pytest.raises(ValueError):
        isotope_mass(isotope)


def test_invalid_number_of_isotopes():
    with pytest.raises(ValueError, match="exactly two"):
        reduced_mass(("H-1",))
