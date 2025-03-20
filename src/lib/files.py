import os


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
    script_dir = os.path.dirname(__file__)
    directory = os.path.join(script_dir, f'../../measurements/{directory}/')
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
    paths = get_measurement_paths(directory)
    return [f'{directory}/{f[:-14]}' for f in paths if f.endswith('-reference.lvm')]