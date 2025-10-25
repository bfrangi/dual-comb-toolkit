from typing import TYPE_CHECKING

import matplotlib.pyplot as plt

if TYPE_CHECKING:
    from typing import Optional

    from numpy.typing import ArrayLike

tight = {"pad": 0.1, "rect": (0.02, 0, 0.99, 0.99)}
article_tight = {"pad": 0.1, "rect": (0.04, 0.07, 0.99, 0.99)}
article_scale = 0.6
article_figsize = (10 * article_scale, 6 * article_scale)


def spectrum_plot(
    wu,
    transmission,
    title,
    xlabel="Waveunit [wu]",
    ylabel="Transmittance [-]",
    **kwargs,
) -> plt:
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
    plt.yscale(kwargs.get("yscale", "linear"))
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
        raise ValueError("Invalid number of arguments")

    color = kwargs.get("color", "b")
    title = kwargs.get("title", "")
    xlabel = kwargs.get("xlabel", "")
    ylabel = kwargs.get("ylabel", "")
    interp = kwargs.get("interp", False)
    size = kwargs.get("size", None)
    save_as = kwargs.get("save_as", None)

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

        filename = f"{os.path.dirname(os.path.realpath(__file__))}/../figures/{save_as}"
        plt.savefig(f"../figures/{save_as}")
        print(f"Figure saved to {filename}")

    return plt


def stem_plot(x: "ArrayLike", y: "ArrayLike", **kwargs) -> plt:
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
    color = kwargs.get("color", "blue")
    title = kwargs.get("title", "")
    xlabel = kwargs.get("xlabel", "")
    ylabel = kwargs.get("ylabel", "")

    from matplotlib import pyplot as plt

    plt.stem(x, y, linefmt=color)
    if title:
        plt.title(title)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)

    return plt


def toggle_series_plot(
    data: "list[tuple[ArrayLike, ArrayLike, str]]",
    title: str,
    xlabel: "Optional[str]" = None,
    ylabel: "Optional[str]" = None,
    cmap=None,
    zoom: float = 1.0,
    figsize: tuple[int] = (9, 6),
    save_path: "Optional[str]" = None,
) -> "tuple[plt.Figure, plt.Axes]":
    """Plot multiple series with a toggleable legend.

    Parameters
    ----------
    data : list of tuples
        Each tuple contains (x, y, label) for the series to plot.
    title : str
        The title of the plot.
    xlabel : str, optional
        The label for the x-axis. Defaults to None.
    ylabel : str, optional
        The label for the y-axis. Defaults to None.
    cmap : matplotlib.colors.Colormap, optional
        A colormap to use for the series. If None, a default color cycle is used.
    zoom : float, optional
        The zoom factor for the figure. Defaults to 1.0.
    figsize : tuple of int, optional
        The size of the figure in inches. Defaults to (9, 6).
    save_path : str, optional
        If provided, the plot will be saved to this path. Defaults to None.

    Returns
    -------
    tuple[plt.Figure, plt.Axes]
        A tuple containing the figure and axes objects.
    """
    fig, ax = plt.subplots(figsize=(figsize[0] / zoom, figsize[1] / zoom), dpi=300)
    series_list = []
    for i, item in enumerate(data):
        x, y, label = item

        if cmap is None:
            (series,) = ax.plot(x, y, "o-", label=label)
        else:
            (series,) = ax.plot(x, y, "o-", label=label, color=cmap(i))

        series_list.append(series)

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout(**tight)

    legend = plt.legend(prop={"size": 6.5})
    legends = legend.get_lines()

    graphs = {}
    for i, leg in enumerate(legends):
        leg.set_picker(True)
        leg.set_pickradius(4)
        graphs[leg] = series_list[i]

    def on_pick(event):
        legend = event.artist
        isVisible = legend.get_visible()

        graphs[legend].set_visible(not isVisible)
        legend.set_visible(not isVisible)

        fig.canvas.draw()

    plt.connect("pick_event", on_pick)

    if save_path:
        plt.savefig(save_path)

    return fig, ax


def config_plot(
    x, y, title, xlabel, ylabel, save_path=None
) -> "tuple[plt.Figure, plt.Axes]":
    """Plot a configuration plot with given x and y data.

    Parameters
    ----------
    x : array-like
        The x values.
    y : array-like
        The y values.
    title : str
        The title of the plot.
    xlabel : str
        The label for the x-axis.
    ylabel : str
        The label for the y-axis.
    save_path : str, optional
        If provided, the plot will be saved to this path. Defaults to None.

    Returns
    -------
    tuple[plt.Figure, plt.Axes]
    fig, ax : tuple
        The figure and axes objects of the plot.
    """
    scale = 1.5
    fig, ax = plt.subplots(figsize=(3 * scale, 2 * scale), dpi=300)
    ax.plot(x, y, "o-", c="b")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.tight_layout(**tight)

    if save_path:
        plt.savefig(save_path)
    return fig, ax


def use_latex():
    """Use LaTeX for rendering text in plots."""
    plt.rcParams.update({"text.usetex": True, "font.family": "Computer Modern"})


cmaps = {
    "CMRmap": 2,
    "brg": 2,
    "winter": 1,
    "twilight": 1,
    "copper": 1,
    "berlin": 1,
    "coolwarm": 1,
    "viridis": 3 / 2,
    "autumn": 3 / 2,
    "summer": 3 / 2,
    "gnuplot": 4 / 3,
}


def get_cmap(name: str, n: int):
    """Get a colormap with n distinct colors.

    Parameters
    ----------
    name : str
        The name of the colormap.
    n : int
        The number of distinct colors needed.

    Returns
    -------
    matplotlib.colors.Colormap
        The requested colormap with n distinct colors.
    """
    if name in cmaps:
        factor = cmaps[name]
        n = int(n * factor)
    return plt.get_cmap(name, n)
