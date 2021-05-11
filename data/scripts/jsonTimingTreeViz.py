#!/usr/bin/env python3

import argparse
import json
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
import sys


def add_pie_piece(ax, fractional_width, fractional_location, color, level, text='', height=0.3, inner_radius=0.5, flatten=False):
    fractional_width *= 2 * np.pi
    fractional_location *= 2 * np.pi
    bar_center = (fractional_location + fractional_width / 2 + 0.05 * fractional_width * level, inner_radius+level*height + height/2)
    annotation_location = (bar_center[0], bar_center[1] + inner_radius)
    if flatten:
        annotation_location = bar_center
    ax.bar(x=fractional_location, width=fractional_width, bottom=inner_radius+level*height, height=height, edgecolor='w', linewidth=1, align='edge', label=text, color=color)
    ax.annotate(text, bar_center, xytext=annotation_location, bbox=dict(boxstyle="round4", fc="w", alpha=0.5),
                arrowprops=dict(color='gray', arrowstyle="-|>", connectionstyle="arc3,rad=-0.2", fc="w"), fontsize='small')

def add_pie_piece_recursively(ax, json_dict, start_location=0., level=0, parent_fractional_width=1., rel_threshold=0.2, abs_threshold=0.05, flatten=False):
    level_sum = sum([json_dict[data]['total'] * 100 for data in json_dict])
    current_location = start_location
    for key in json_dict:
        data = json_dict[key]
        relative_fractional_width = data['total'] * 100 / level_sum
        fractional_width = relative_fractional_width * parent_fractional_width
        if relative_fractional_width >= rel_threshold and fractional_width >= abs_threshold:
            legend_text = key + '\n' + "rel: {0:.2f}".format(round(relative_fractional_width * 100, 2)) + '%\n' + "abs: {0:.2f}".format(round(fractional_width * 100, 2)) + '%\n' + \
                          "total time: {0:.2f}".format(data['total'])
            add_pie_piece(ax, fractional_width, current_location, cm.Reds(fractional_width), level, text=legend_text, flatten=flatten)
            add_pie_piece_recursively(ax, data['childs'], start_location=current_location, level=level+1, parent_fractional_width=fractional_width, rel_threshold=rel_threshold, abs_threshold=abs_threshold, flatten=flatten)
        current_location += fractional_width

def main():
    parser = argparse.ArgumentParser(description='Tool to plot walberla TimingTree JSON output files as pie chart.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('file', type=str,
                        help='JSON file generated from a walberla TimingTree')
    parser.add_argument('-r', '--rel-threshold', type=float, default=0.2,
                        help='threshold for the relative amount of time spent to be included in the plot')
    parser.add_argument('-a', '--abs-threshold', type=float, default=0.05,
                        help='threshold for the absolute amount of time spent to be included in the plot')
    parser.add_argument('-o', '--output-file', type=str,
                        help='writes the pie chart to file')
    parser.add_argument('--legend', action='store_true',
                        help='plot traditional color legend')
    parser.add_argument('--show-filename', action='store_true',
                        help='displays the passed file name in the image')
    parser.add_argument('--flatten', action='store_true',
                        help='flattening out pie')

    args = parser.parse_args()

    if not os.path.exists(args.file):
        sys.exit('Input file does not exist.')

    fig, ax = plt.subplots(subplot_kw=dict(polar=(not args.flatten)))

    with open(args.file) as f:
        data = json.load(f)

        add_pie_piece_recursively(ax, data, rel_threshold=args.rel_threshold, abs_threshold=args.abs_threshold, flatten=args.flatten)

        if args.show_filename:
            plt.title(args.file, y=1.0)

        ax.set_axis_off()

        if args.legend:
            plt.legend(loc='best')

        if args.output_file:
            plt.savefig(args.output_file)
        else:
            plt.show()

if __name__ == '__main__':
    main()

