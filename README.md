# Dual Comb Processing

Useful scripts for processing dual comb data.

## Setup

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
chmod +x ./radis-patch/patch.sh
./radis-patch/patch.sh # Optional, to suppress vaex output
```

4. [Optional] If you want to use LaTeX for plotting, install the `texlive` package:

```bash
sudo apt-get update
sudo apt-get install texlive-latex-extra texlive-fonts-extra dvipng
```

## Using GPU Acceleration

Check out the `DEVICE_ID` setting in `src/lib/gpu.py`. By default, it is set to `"nvidia"` which will
use the NVIDIA GPU if available. If you want to use a different GPU, run the `identify-gpu.py` script to
list available GPUs and set the `DEVICE_ID` accordingly (it can be the number of the GPU in the 
output list or a string contained in the name of the device).