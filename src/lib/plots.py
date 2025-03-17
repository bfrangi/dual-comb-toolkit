import matplotlib.pyplot as plt

tight = {
    'pad': 0.1,
    'rect': (-0.02, 0, 0.99, 0.99)
}


def spectrum_plot(wu, transmission, title, xlabel='Waveunit [wu]', ylabel='Transmittance [-]') -> plt:
    """
    Plot the given transmission spectrum.

    Parameters
    ----------
    wu : np.ndarray
        The waveunit array (in Hz, cm⁻¹, nm, ...).
    transmission : np.ndarray
        The transmission spectrum.
    title : str
        The title of the plot.
    xlabel : str
        The x-axis label.
    ylabel : str
        The y-axis label.

    Returns
    -------
    plt
        The matplotlib plot of the transmission spectrum.
    """
    plt.plot(wu, transmission)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout(**tight)
    return plt

