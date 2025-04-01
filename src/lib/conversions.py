from lib.constants import c


# Pressure conversions


def pa_to_bar(p: float) -> float:
    """
    Convert pressure from Pascals to bars.

    Parameters
    ----------
    p : float
        Pressure in Pa.

    Returns
    -------
    float
        Pressure in bar.
    """
    from astropy import units as u
    return (p * u.Pa).to(u.bar).value


# Frequency, wavelength and wavenumber conversions


def wavelength_to_wavenumber(wl: float) -> float:
    """
    Convert wavelength to wavenumber.

    Parameters
    ----------
    wl : float
        Wavelength in nm.

    Returns
    -------
    float
        Wavenumber in cm^-1.
    """
    return 1 / wl * 1e7


def wavenumber_to_wavelength(wn: float) -> float:
    """
    Convert wavenumber to wavelength.

    Parameters
    ----------
    wn : float
        Wavenumber in cm^-1.

    Returns
    -------
    float
        Wavelength in nm.
    """
    return 1 / wn * 1e7


def delta_wavelength_to_delta_wavenumber(delta_wl, central_wl):
    """
    Convert a small wavelength range around a central wavelength to a wavenumber 
    range. This is done by using the formula: 
    `delta_wn = delta_wl / central_wl^2`.

    Parameters
    ----------
    delta_wl : float
        The wavelength range in nm.
    central_wl : float
        The central wavelength in nm.

    Returns
    -------
    float
        The wavenumber range in cm^-1.
    """
    return (delta_wl * 1e-9) / (central_wl * 1e-9) ** 2 * 1e-2


def delta_frequency_to_delta_wavelength(delta_fr, central_fr):
    """
    Convert a small frequency range around a central frequency to a wavelength
    range. This is done by using the formula:
    `delta_wl = delta_fr * c / central_fr^2`. 

    Parameters
    ----------
    delta_fr : float
        The frequency range in Hz.
    central_wl : float
        The central frequency in Hz.

    Returns
    -------
    float
        The wavelength range in nm.
    """
    return delta_fr * c / central_fr ** 2 * 1e9
