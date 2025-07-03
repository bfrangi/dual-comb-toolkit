#!/usr/bin/env python3

import platform
import shutil
import subprocess
import sys

this_file_dir = (
    __file__.rsplit("/", 1)[0] if "/" in __file__ else __file__.rsplit("\\", 1)[0]
)
python = sys.executable
script = f"{this_file_dir}/fit-simulated-measurement.py"


def launch_command_windows(base_config, n):
    for i in range(n):
        cmd = f'"{python}" "{script}" --config {base_config}_{i}.txt'
        subprocess.Popen(["start", "cmd", "/k", cmd], shell=True)
        print(f"Launched terminal {i} with command: {cmd}")


def launch_command_linux(base_config, n):
    terminal = None
    terminals = ["gnome-terminal", "xterm", "konsole", "xfce4-terminal", "lxterminal"]
    for t in terminals:
        if shutil.which(t):
            terminal = t
            break

    print(f"No supported terminal found. Install one of: {', '.join(terminals)}.")
    if not terminal:
        print(f"No supported terminal found. Install one of: {', '.join(terminals)}.")
        sys.exit(1)

    for i in range(n):
        cmd = f'"{python}" "{script}" --config {base_config}_{i}.txt'
        if terminal == "gnome-terminal":
            subprocess.Popen([terminal, "--", "bash", "-c", f"{cmd}; exec bash"])
        elif terminal == "xterm":
            subprocess.Popen([terminal, "-e", f"{cmd}; bash"])
        elif terminal == "konsole":
            subprocess.Popen([terminal, "-e", cmd])
        elif terminal == "xfce4-terminal":
            subprocess.Popen([terminal, "-e", f"bash -c '{cmd}; exec bash'"])
        elif terminal == "lxterminal":
            subprocess.Popen([terminal, "-e", cmd])
        print(f"Launched terminal {i} with command: {cmd}")


def main():
    if len(sys.argv) != 3:
        print(
            "Usage: python launch_simulated_fittings.py base_configurations_name number_of_configurations"
        )
        sys.exit(1)

    base_config = sys.argv[1]
    n = int(sys.argv[2])
    system = platform.system()

    if system == "Windows":
        launch_command_windows(base_config, n)
    elif system == "Linux":
        launch_command_linux(base_config, n)
    else:
        print(f"Unsupported OS: {system}")
        sys.exit(1)


if __name__ == "__main__":
    main()
