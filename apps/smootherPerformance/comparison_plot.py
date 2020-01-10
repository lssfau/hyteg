import numpy as np
import argparse
import csv
from matplotlib import pyplot as plt


def read_col(filename, colname):
    l2_error_list = []
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            l2_error_list.append(float(row[colname]))
    return l2_error_list


def read_res(filename):
    return read_col(filename, 'resNorm')


parser = argparse.ArgumentParser(description='Plot the residual rates of different simulations.')
parser.add_argument(
    '--file-paths',
    help='the filepaths to the csv files',
    nargs='+',
    type=str,
    required=True
)
parser.add_argument(
    '--labels',
    help='the labels which appear in the plot',
    nargs='+',
    type=str,
    required=False
)
parser.add_argument(
    '--title',
    help='the title of the plot',
    type=str,
    default=''
)
parser.add_argument(
    '--max-iter',
    help='the last iteration we should take into account',
    type=int,
    required=False
)
parser.add_argument(
    '--y-limit',
    help='should the y axis be restricted?',
    type=float,
    required=False
)
PLOT_RESIDUAL = 'residual'
PLOT_RESIDUAL_RATE = 'residual-rate'
parser.add_argument('--plot', help='What should we plot?', type=str, default=PLOT_RESIDUAL,
                    choices=[PLOT_RESIDUAL, PLOT_RESIDUAL_RATE])

pargs = parser.parse_args()

l2_all_errors_list = []
for file_path in pargs.file_paths:
    if pargs.plot == PLOT_RESIDUAL:
        entries = read_res(file_path)
    elif pargs.plot == PLOT_RESIDUAL_RATE:
        entries = np.array(read_res(file_path))
        entries = entries[1:] / entries[0:-1]
    else:
        print('plot command not recognized')
        exit()
    l2_all_errors_list.append(entries)

iteration_list = list(range(1, len(l2_all_errors_list[0]) + 1))

for idx, (label, l2_error_list) in enumerate(zip(pargs.labels, l2_all_errors_list)):
    sym = ('^:', '.:', '+:', 'x:')[idx % 4]
    if pargs.plot == PLOT_RESIDUAL_RATE:
        plt.plot(iteration_list[0:pargs.max_iter], l2_error_list[0:pargs.max_iter], sym, label=label)
    else:
        plt.semilogy(iteration_list[0:pargs.max_iter], l2_error_list[0:pargs.max_iter], sym, label=label)

plt.xlabel('Multigrid iterations')
if pargs.plot == PLOT_RESIDUAL:
    plt.ylabel('residual $r$')
elif pargs.plot == PLOT_RESIDUAL_RATE:
    plt.ylabel('residual rate $\\frac{r_i}{r_{i-1}}$')
else:
    print('plot command not recognized')
    exit()

if pargs.y_limit:
    plt.ylim([0, pargs.y_limit])

plt.legend()
plt.title(pargs.title)
plt.grid(True)
plt.show()
