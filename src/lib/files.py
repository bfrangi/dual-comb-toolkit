import csv
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
    return os.path.normpath(os.path.join(script_dir, "../../"))


def get_figures_path() -> str:
    """
    Get the figures folder path.

    Returns
    -------
    str
        The figures folder path.
    """
    root_path = get_root_path()
    return os.path.join(root_path, "figures/")


def get_reports_path() -> str:
    """
    Get the reports folder path.

    Returns
    -------
    str
        The reports folder path.
    """
    root_path = get_root_path()
    return os.path.join(root_path, "reports/")


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
    directory = os.path.join(root_path, f"measurements/{directory}/")
    return [f for f in os.listdir(directory) if f.endswith(".lvm")]


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
    return [
        f"{directory}/{f[:-14]}"
        for f in measurement_paths
        if f.endswith("-reference.lvm")
    ]


def initialize_csv_report(filename: str, headers: tuple[str]) -> None:
    """
    Initialize a CSV report file with the given headers.

    Parameters
    ----------
    filename : str
        The name of the CSV file to create.
    headers : list[str]
        The headers for the CSV file.
    """
    csv_path = f"{get_reports_path()}{filename}.csv"

    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile, delimiter=",")
        writer.writerow(headers)


def append_to_csv_report(filename: str, data: tuple) -> None:
    csv_path = f"{get_reports_path()}{filename}.csv"

    with open(csv_path, "a", newline="") as csvfile:
        writer = csv.writer(csvfile, delimiter=",")
        writer.writerow(data)


def create_figures_folder(foldername) -> str:
    """
    Create a folder in the figures directory if it does not exist.

    Parameters
    ----------
    foldername : str
        The name of the folder to create in the figures directory.

    Returns
    -------
    str
        The path to the created folder.
    """
    folder_path = os.path.join(get_figures_path(), foldername)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    return folder_path
