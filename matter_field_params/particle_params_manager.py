from .scipy_codata import get_codata_table


SUPPORTED_PARTICLES = ("electron", "proton", "deuteron", "muon", "positron")

PARTICLE_CHARGES = {
    "electron": -1.0,
    "proton": 1.0,
    "deuteron": 1.0,
    "muon": -1.0,
    "positron": 1.0,
}

PREDEFINED_3BODY_SYSTEMS = {
    # Tuple layout: (reference_particle, particle1, particle2)
    "Ps^-": ("electron", "positron", "electron"),
    "H2+": ("proton", "proton", "electron"),
    "HD+": ("deuteron", "proton", "electron"),
    "D2+": ("deuteron", "deuteron", "electron"),
}

FRAME_ALIASES = {
    "reference_center": "reference_center",
    "reference_frame": "reference_center",
    "rc": "reference_center",
    "nuclear_center": "nuclear_center",
    "nc": "nuclear_center",
    "ncm": "nuclear_center",
}


def _normalize_frame(frame: str) -> str:
    try:
        return FRAME_ALIASES[frame.lower()]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported frame: '{frame}'. Choose between 'reference_center' "
            "(aliases: 'reference_frame', 'rc') and 'nuclear_center' "
            "(aliases: 'nc', 'ncm')."
        ) from exc


def predefined_3body_systems(system: str) -> tuple[str, str, str]:
    """Return a system tuple as (reference_particle, particle1, particle2)."""
    aliases = {
        "ps-": "Ps^-",
        "ps^-": "Ps^-",
        "h2+": "H2+",
        "h2p": "H2+",
        "hd+": "HD+",
        "hdp": "HD+",
        "d2+": "D2+",
        "d2p": "D2+",
    }
    normalized = aliases.get(system.lower(), system)

    try:
        return PREDEFINED_3BODY_SYSTEMS[normalized]
    except KeyError as exc:
        supported = ", ".join(PREDEFINED_3BODY_SYSTEMS.keys())
        raise ValueError(f"Unsupported system '{system}'. Choose from: {supported}") from exc


def get_reduced_masses_and_charges_3body(
    particles=("proton", "proton", "electron"),
    codata_set="current",
    frame="nuclear_center",
):
    """Compute reduced masses and reduced charges for a generic 3-body system.

    Parameters
    ----------
    frame : str
        One of "reference_center" or "nuclear_center".
        Aliases: "reference_frame"/"rc" and "nc"/"ncm".
    particles : tuple[str, str, str]
        Tuple on the form (reference_particle, particle1, particle2).
    codata_set : str
        The CODATA set to use for physical constants. Can be "current" or one
        of the supported CODATA years: "2006", "2010", "2014", "2018", "2022".

    Notes
    -----
    Supported particles: electron, proton, deuteron, muon, positron.
    Masses are in atomic units (electron mass units).
    """

    try:
        reference_particle, particle1, particle2 = particles
    except Exception as exc:
        raise ValueError(
            "particles must be a tuple/list of exactly 3 elements: "
            "(reference_particle, particle1, particle2)."
        ) from exc

    constants_table, codata_version = get_codata_table(codata_set)
    frame = _normalize_frame(frame)

    for role, particle in (
        ("reference_particle", reference_particle),
        ("particle1", particle1),
        ("particle2", particle2),
    ):
        if particle not in SUPPORTED_PARTICLES:
            raise ValueError(
                f"Unsupported {role}: '{particle}'. Choose from {SUPPORTED_PARTICLES}."
            )

    def _mass_au(particle: str) -> float:
        if particle in ("electron", "positron"):
            return 1.0

        key = f"{particle}-electron mass ratio"
        try:
            return constants_table[key][0]
        except KeyError as exc:
            raise ValueError(
                f"Mass ratio '{key}' not available in CODATA {codata_version}."
            ) from exc

    m_ref = _mass_au(reference_particle)
    m_1 = _mass_au(particle1)
    m_2 = _mass_au(particle2)

    m_tot = m_ref + m_1 + m_2

    q_ref = PARTICLE_CHARGES[reference_particle]
    q_1 = PARTICLE_CHARGES[particle1]
    q_2 = PARTICLE_CHARGES[particle2]

    q_tot = q_ref + q_1 + q_2

    mu_1 = (m_ref * m_1) / (m_ref + m_1)

    if frame == "reference_center":
        mu_2 = (m_ref * m_2) / (m_ref + m_2)
        xi_1 = q_1 - (m_1 / m_tot) * q_tot
    else:
        m_pair = m_ref + m_1
        q_pair = q_ref + q_1
        mu_2 = (m_pair * m_2) / (m_pair + m_2)
        xi_1 = q_1 - (m_1 / m_pair) * q_pair

    xi_2 = q_2 - (m_2 / m_tot) * q_tot

    return mu_1, mu_2, xi_1, xi_2, m_ref, m_1, m_2, q_ref, q_1, q_2
