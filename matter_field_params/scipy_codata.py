from scipy.constants import physical_constants
from scipy.constants import _codata as scipy_codata


SUPPORTED_CODATA_SETS = ("2006", "2010", "2014", "2018", "2022")


def _resolve_current_codata_year() -> str:
    reference_keys = ("electron mass", "proton mass", "deuteron mass")

    for year in reversed(SUPPORTED_CODATA_SETS):
        table = getattr(scipy_codata, f"_physical_constants_{year}")
        if all(physical_constants[key][0] == table[key][0] for key in reference_keys):
            return year

    return "current"

def get_codata_table(codata_set: str):
    if codata_set == "current":
        return physical_constants, _resolve_current_codata_year()

    if codata_set not in SUPPORTED_CODATA_SETS:
        supported = ", ".join(["current", *SUPPORTED_CODATA_SETS])
        raise ValueError(f"Unsupported codata_set '{codata_set}'. Use one of: {supported}")

    return getattr(scipy_codata, f"_physical_constants_{codata_set}"), codata_set


if __name__ == "__main__":
    previous_set = None
    previous_proton_mass_au = None
    previous_deuteron_mass_au = None

    for codata_set in SUPPORTED_CODATA_SETS:
        table, resolved_set = get_codata_table(codata_set)

        electron_mass = table["electron mass"][0]
        proton_mass_au = table["proton mass"][0] / electron_mass
        deuteron_mass_au = table["deuteron mass"][0] / electron_mass

        print(f"CODATA set      : {resolved_set}")
        print(f"M(proton)   (au): {proton_mass_au:.8f}")
        print(f"M(deuteron) (au): {deuteron_mass_au:.8f}")

        if previous_set is not None:
            delta_proton = proton_mass_au - previous_proton_mass_au
            delta_deuteron = deuteron_mass_au - previous_deuteron_mass_au
            print(f"Delta M(proton)   vs {previous_set}: {delta_proton:+.2g} au")
            print(f"Delta M(deuteron) vs {previous_set}: {delta_deuteron:+.2g} au")

        previous_set = resolved_set
        previous_proton_mass_au = proton_mass_au
        previous_deuteron_mass_au = deuteron_mass_au
        print()