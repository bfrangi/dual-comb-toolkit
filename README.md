# Dual-Comb Toolkit

Useful scripts for working with dual-comb spectroscopy (DCS), including:
- Fast simulation of absorption spectra using the [Radis](https://github.com/radis/radis) library,
  compatible with GPU.
- Processing of dual-comb sample/reference measurements (baseline correction, normalization, 
  RF-to-optical mapping, noisy tooth filtering, etc).
- Fitting of processed dual-comb spectra to retrieve gas concentration.
- Mapping of concentration across multiple measurement positions to create 2D maps.
- Analysis of temporal concentration evolution within measurements. 
- Simulation of real, noisy measurements, for optimization of experimental configurations.

## Setup

Before starting, if you are going to be using the HITRAN/HITEMP databases, create a Hitran
login [here](https://hitran.org/login/). Other databases may require their own credentials.

Once you have your credentials, you can proceed with the setup below, depending on your
operating system.

<details>
<summary><b>Linux</b></summary>

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
   ```

4. [Optional] If you want to use LaTeX for plotting, install the `texlive` package:

   ```bash
   sudo apt-get update
   sudo apt-get install texlive-latex-extra texlive-fonts-extra dvipng
   ```

</details>

<details>
<summary><b>Windows</b></summary>

1. Clone the repository:

   ```bash
   git clone https://github.com/bfrangi/dual-comb-toolkit.git
   cd dual-comb-toolkit
   ```

2. Install `python3.12` from the Microsoft Store.

3. Create a virtual environment:

   ```bash
   python3.12 -m venv .venv
   .\.venv\Scripts\activate.bat
   ```

4. Install dependencies and apply the radis patch (this second part is optional):

   ```bash
   pip install -r requirements-win.txt
   ```

5. [Optional] If you want to use LaTeX for plotting, install the `texlive` package by downloading
   the installer from [here](https://www.tug.org/texlive/windows.html) or the ISO from 
   [here](https://www.tug.org/texlive/acquire-iso.html) and following the installation
   instructions.

</details>

## Usage

If using HITRAN/HITEMP, you will need to log in with your HITRAN credentials the first time you
simulate an absorption spectrum. Other databases may require their own credentials.

**Note**: The first time you run the simulator, the database will be downloaded, and this could
take some time. Please be patient!

<details>
<summary><b>Line simulation</b></summary>

You can simulate absorption spectra using the `src/line-simulator.py` script. You can define the
simulation parameters directly in the script. Navigate to the `src` directory and run the script
as follows:

```bash
python line-simulator.py
```

</details>

<details>
<summary><b>Manual measurement processing</b></summary>

You can manually process dual-comb measurements using the `src/process-measurement.py` script. You
can define the measurement parameters directly in the script, including the path to the measurement 
to process and the path to the set of measurements to use for baseline correction. This last part is
useful to remove etalon effects and other systematic noise from the measurement. Navigate to the 
`src` directory and run the script as follows:

```bash
python process-measurement.py
```

Here is an example of the output plot:

![Example processed plot](assets/process-measurement-example.svg)

To characterize the baseline, there is also the `src/process-baseline.py` script. You can run it as
follows:

```bash
python process-baseline.py
```

And there is also a simple script to view the raw spectrum of any measurement, `src/measurement-spectrum.py`:

```bash
python measurement-spectrum.py
```

An example of the output plot from this script is shown below:

![Spectrum](assets/measurement-spectrum-example.svg)

</details>

<details>
<summary><b>Concentration fitting</b></summary>

You can fit processed dual-comb measurements to retrieve gas concentration using the
`src/fit-measurement.py` script. You can define the fitting parameters directly in the script,
very similarly to the `src/process-measurement.py` script. Navigate to the `src` directory and run
the script as follows:

```bash
python fit-measurement.py
```

Here is an example of the output plot:

![Example fit plot](assets/fit-measurement-example.svg)

</details>

<details>
<summary><b>Concentration fitting</b></summary>

</details>

## Using GPU Acceleration

Check out the `GPU_DEVICE_ID` setting in `src/lib/defaults.py`. By default, it is set to `"nvidia"`
which will use the NVIDIA GPU if available. If you want to use a different GPU, run the
`src/identify-gpu.py` script to list available GPUs and set the `GPU_DEVICE_ID` accordingly (it can 
be the number of the GPU in the output list or a string contained in the name of the device).

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

## Acknowledgements

This project uses the [RADIS](https://github.com/radis/radis) library for
high-resolution molecular spectroscopy calculations.
