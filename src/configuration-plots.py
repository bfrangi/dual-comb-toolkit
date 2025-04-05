import csv

import matplotlib.cm as cm
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure

from lib.files import get_figures_path, get_reports_path
from lib.plots import tight

report_name = 'report'

report_path = f"{get_reports_path()}{report_name}.csv"

with open(report_path) as csvfile:
    reader = csv.reader(csvfile)
    rows = list(reader)[1:]

data: dict[int, dict[str, list[float]]] = {}

for row in rows:
    n_teeth = int(row[0])
    spacing = float(row[1])
    concentration = float(row[2])
    sdv = float(row[3])

    if n_teeth not in data:
        data[n_teeth] = {
            'spacings': [],
            'concentrations': [],
            'sdvs': []
        }

    data[n_teeth]['spacings'].append(spacing/1e9)
    data[n_teeth]['concentrations'].append(concentration)
    data[n_teeth]['sdvs'].append(sdv)


for n_teeth, d in data.items():
    plt.plot(d['spacings'], d['sdvs'], 'o-', c='b')
    plt.title('Standard deviation of the concentration as a function\n' +
              f'of the comb spacing for {n_teeth} teeth')
    plt.xlabel('Spacing (GHz)')
    plt.ylabel('Standard deviation (VMR)')
    plt.tight_layout(**tight)
    plt.savefig(f'{get_figures_path()}conf-sdv-{n_teeth}.svg')
    plt.clf()

    plt.plot(d['spacings'], d['concentrations'], 'o-', c='b')
    plt.title(f'Concentration as a function of the comb spacing for {n_teeth} teeth')
    plt.xlabel('Spacing (GHz)')
    plt.ylabel('Concentration (VMR)')
    plt.tight_layout(**tight)
    plt.savefig(f'{get_figures_path()}conf-conc-{n_teeth}.svg')
    plt.clf()

num_series = len(data)
cmap = cm.get_cmap('viridis', num_series)

figure(figsize=(9, 6), dpi=80)

for i, item in enumerate(data.items()):
    n_teeth = item[0]
    d = item[1]
    plt.plot(d['spacings'], d['sdvs'], 'o-', label=f'{n_teeth} teeth', color=cmap(i))

plt.title('Standard deviation of the concentration as a function\n' +
          'of the comb configuration')
plt.xlabel('Spacing (GHz)')
plt.ylabel('Standard deviation (VMR)')
plt.legend()
plt.tight_layout(**tight)
plt.savefig(f'{get_figures_path()}conf-sdv.svg')
plt.clf()

for i, item in enumerate(data.items()):
    n_teeth = item[0]
    d = item[1]
    plt.plot(d['spacings'], d['concentrations'], 'o-', label=f'{n_teeth} teeth', color=cmap(i))

plt.title('Concentration as a function of the comb configuration')
plt.xlabel('Spacing (GHz)')
plt.ylabel('Concentration (VMR)')
plt.legend()
plt.tight_layout(**tight)
plt.savefig(f'{get_figures_path()}conf-conc.svg')
plt.clf()


