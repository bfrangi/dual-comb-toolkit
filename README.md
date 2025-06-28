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
python3 -m venv .venv
source .venv/bin/activate
```

3. Install dependencies:

```bash
pip install -r requirements.txt
```

4. [Optional] If you want to use LaTeX for plotting, install the `texlive` package:

```bash
sudo apt-get update
sudo apt-get install texlive-latex-extra texlive-fonts-extra
```