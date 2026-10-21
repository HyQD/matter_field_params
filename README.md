# matter_field_params

Utilities for matter and electromagnetic-field parameters.

## Isotope masses

Isotopes must be specified explicitly. Masses can be returned in unified
atomic mass units or electron-mass atomic units:

```python
from matter_field_params import isotope_mass, nuclear_mass, reduced_mass

atomic_li7 = isotope_mass("Li-7", units="u")
nuclear_li7 = nuclear_mass("Li-7", units="electron_mass")
mu_lih = reduced_mass(("Li-7", "H-1"), units="electron_mass")
```

`isotope_mass` returns the tabulated neutral-atom mass by default. The
`nuclear_mass` function and the default behavior of `reduced_mass` subtract
the electron rest masses to approximate bare-nucleus masses. Electronic
binding energies are not included in this approximation. Use
`mass_kind="atomic"` with `reduced_mass` when neutral-atom masses are desired.
