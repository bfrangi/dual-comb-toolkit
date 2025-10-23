# Dual Comb Toolkit

Useful scripts for processing dual comb data.

## Setup

### Linux

1. Clone the repository:

```bash
git clone https://github.com/bfrangi/dual-comb-toolkit.git
cd dual-comb-toolkit
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
git clone https://github.com/bfrangi/dual-comb-toolkit.git
cd dual-comb-toolkit
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

## Issues with LaTeX (on Linux)

If you encounter `FileNotFoundError: Matplotlib's TeX implementation searched for a file named 'cmr10.tfm' in your texmf tree, but could not find it.`, just run the [following command](https://stackoverflow.com/a/79243265/15159198):

```bash
sudo mv /usr/bin/luatex /usr/bin/luatex.bk 
```

If you encounter `LaTeX Error: File 'type1ec.sty' not found.`, install the `cm-super` package (see [here](https://github.com/matplotlib/matplotlib/issues/16911)):

```bash
sudo apt install cm-super
```

If you encounter `RuntimeError: Failed to process string with tex because dvipng could not be found`, install the `dvipng` package:

```bash
sudo apt install dvipng
```