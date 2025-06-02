from typing import TYPE_CHECKING

import matplotlib.pyplot as plt

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

tight = {
    'pad': 0.1,
    'rect': (0.02, 0, 0.99, 0.99)
}


def spectrum_plot(wu, transmission, title, xlabel='Waveunit [wu]', ylabel='Transmittance [-]', **kwargs) -> plt:
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

    Other Parameters
    ----------------
    yscale : str, optional
        The scale of the y-axis. Default is 'linear'.

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
    plt.yscale(kwargs.get('yscale', 'linear'))
    return plt


def scatter_plot(*args, **kwargs) -> plt:
    """Plot a scatter plot of x and y with optional interpolation.

    Parameters
    ----------
    x : array, optional
        The x values.
    y : array
        The y values.

    Keyword Arguments
    -----------------
    color : str
        The color of the plot.
    title : str
        The title of the plot.
    xlabel : str
        The x-axis label.
    ylabel : str
        The y-axis label.
    interp : bool
        Whether to interpolate the data.
    size : tuple
        The size of the plot.
    save_as : str
        The path to save the figure.

    Returns
    -------
    matplotlib.pyplot
        The generated plot.
    """
    if len(args) == 1:
        x = range(len(args[0]))
        y = args[0]
    elif len(args) == 2:
        x, y = args
    else:
        raise ValueError('Invalid number of arguments')

    color = kwargs.get('color', 'b')
    title = kwargs.get('title', '')
    xlabel = kwargs.get('xlabel', '')
    ylabel = kwargs.get('ylabel', '')
    interp = kwargs.get('interp', False)
    size = kwargs.get('size', None)
    save_as = kwargs.get('save_as', None)

    if type(size) not in [list, tuple] or len(size) != 2:
        size = None

    plt.scatter(x, y, c=color)

    if interp:
        from lib.fitting import loose_interpolation
        x, y = loose_interpolation(x, y)
        plt.plot(x, y, c=color)
    if title:
        plt.title(title)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if size:
        plt.gcf().set_size_inches(size[0], size[1])
    if save_as:
        import os
        filename = f'{os.path.dirname(os.path.realpath(__file__))}/../figures/{save_as}'
        plt.savefig(f'../figures/{save_as}')
        print(f'Figure saved to {filename}')

    return plt


def stem_plot(x: 'ArrayLike', y: 'ArrayLike', **kwargs) -> plt:
    """Plot a stem plot of x and y.

    Parameters
    ----------
    x : array-like
        The x values.
    y : array-like
        The y values.

    Keyword Arguments
    -----------------
    color : str, optional
        The color of the plot.
    title : str, optional
        The title of the plot.
    xlabel : str, optional
        The x-axis label.
    ylabel : str, optional
        The y-axis label.

    Returns
    -------
    matplotlib.pyplot
        The generated plot.
    """
    color = kwargs.get('color', 'blue')
    title = kwargs.get('title', '')
    xlabel = kwargs.get('xlabel', '')
    ylabel = kwargs.get('ylabel', '')

    from matplotlib import pyplot as plt

    plt.stem(x, y, linefmt=color)
    if title:
        plt.title(title)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)

    return plt
