import csv
import itertools
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import pylab
import sys

def read_csv(csv_file):
    data = []
    with open(csv_file, mode='r') as infile:
        reader = csv.DictReader(infile)
        for i, r in enumerate(reader):
            if i == 0:
                info = r
            else:
                data.append(r)
    return info, data

def plot_comparison(csv_file, cell_id, outer_iteration):
    info, data = read_csv(csv_file)
    num_processes = info['num_processes']
    plot_data = [float(v["time"]) for v in data if v["outer_iteration"] == str(outer_iteration) and v["cell_id"] == str(cell_id)]
    # plt.subplot(1, len(db_files), i+1)
    # plt.plot(plot_data, label="timings");
    # plt.plot(sorted(plot_data), label="timings sorted");
    plt.bar(list(range(len(plot_data))), plot_data, width=1.0, align='center')
    plt.title("{} processes, cell {}, outer iteration {}".format(int(num_processes), cell_id, outer_iteration))
    plt.xlabel("ranks")
    plt.ylabel("time in seconds")
    plt.show()

if __name__ == "__main__":
    usage = "Plot snoop filter benchmark csv data.\nUsage: plot_csv.py <data.csv> <cell_id> <outer_iteration>"
    if len(sys.argv) != 4:
        print(usage)
    else:
        plot_comparison(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))

    