# Dual Comb Toolbox

Useful scripts for processing dual comb data.

## Setup

### Linux

1. Clone the repository:

```bash
git clone https://github.com/bfrangi/dual-comb-processing.git
cd dual-comb-processing
```

2. Install `virtualenv` and create a virtual environment:

```bash
sudo apt-get update
sudo apt-get install python3-venv
python3.10 -m venv .venv
source .venv/bin/activate
```

3. Install dependencies and apply the radis patch (this second part is optional):

```bash
pip install -r requirements.txt
python radis-patch/apply.py # Optional, to suppress vaex output
```

4. [Optional] If you want to use LaTeX for plotting, install the `texlive` package:

```bash
sudo apt-get update
sudo apt-get install texlive-latex-extra texlive-fonts-extra dvipng
```

### Windows

1. Clone the repository:

```bash
git clone https://github.com/bfrangi/dual-comb-processing.git
cd dual-comb-processing
```

2. Install `Python3.10` from the Microsoft Store.

3. Create a virtual environment:

```bash
python3.10 -m venv .venv
.\.venv\Scripts\activate.bat
```

4. Install dependencies and apply the radis patch (this second part is optional):

```bash
pip install -r requirements-win.txt
python radis-patch\apply.py # Optional, to suppress vaex output
```

5. [Optional] If you want to use LaTeX for plotting, install the `texlive` package by downloading
   the installer from [here](https://www.tug.org/texlive/windows.html) or the ISO from 
   [here](https://www.tug.org/texlive/acquire-iso.html) and following the installation
   instructions.

## Using GPU Acceleration

Check out the `DEVICE_ID` setting in `src/lib/gpu.py`. By default, it is set to `"nvidia"` which will
use the NVIDIA GPU if available. If you want to use a different GPU, run the `identify-gpu.py` script to
list available GPUs and set the `DEVICE_ID` accordingly (it can be the number of the GPU in the 
output list or a string contained in the name of the device).

## Issues with plots

If you encounter `UserWarning: FigureCanvasAgg is non-interactive, and thus cannot be shown plt.show()`,
just install `PyQt6`:

```bash
pip install PyQt6
```
