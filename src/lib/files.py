import os


def get_root_path() -> str:
    """
    Get the root project path.

    Returns
    -------
    str
        The root project path.
    """
    script_dir = os.path.dirname(__file__)
    return os.path.normpath(os.path.join(script_dir, '../../'))


def get_figures_path() -> str:
    """
    Get the figures folder path.

    Returns
    -------
    str
        The figures folder path.
    """
    root_path = get_root_path()
    return os.path.join(root_path, 'figures/')


def get_reports_path() -> str:
    """
    Get the reports folder path.

    Returns
    -------
    str
        The reports folder path.
    """
    root_path = get_root_path()
    return os.path.join(root_path, 'reports/')


def get_measurement_paths(directory: str) -> list[str]:
    """
    Get the measurement paths from the directory.

    Parameters
    ----------
    directory : str
        The directory containing the measurements.

    Returns
    -------
    list[str]
        The measurement paths.
    """
    root_path = get_root_path()
    directory = os.path.join(root_path, f'measurements/{directory}/')
    return [f for f in os.listdir(directory) if f.endswith('.lvm')]


def get_measurement_names(directory: str) -> list[str]:
    """
    Get the measurement names from the directory.

    Parameters
    ----------
    directory : str
        The directory containing the measurements.

    Returns
    -------
    list[str]
        The names of the measurements.
    """
    measurement_paths = get_measurement_paths(directory)
    return [f'{directory}/{f[:-14]}' for f in measurement_paths if f.endswith('-reference.lvm')]
