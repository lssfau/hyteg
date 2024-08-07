#!/usr/bin/env python3

import argparse
import json

import numpy as np
import pandas as pd
import os
import sys
from pathlib import Path


def pad_string_with_spaces(string: str, width: int):
    padding = ' ' * (width - len(string))
    return string + padding


def replace_whitespaces(string: str, old: str = ' ', new: str = '-'):
    return string.replace(old, new)


def fill_table_recursively(json_dict: dict, lists: dict,
                           start_location: float = 0., level: int = 0, parent_fractional_width: float = 1.,
                           parent_keys: str = '',
                           rel_threshold: float = 0.2, abs_threshold: float = 0.05, max_depth: int = float('inf')):
    if level <= max_depth:
        level_sum = sum([json_dict[data]['total'] * 100 for data in json_dict])
        current_location = start_location
        for key in json_dict:
            data = json_dict[key]
            current_key_chain = parent_keys + '.' + key
            relative_fractional_width = data['total'] * 100 / level_sum
            fractional_width = relative_fractional_width * parent_fractional_width
            if relative_fractional_width >= rel_threshold and fractional_width >= abs_threshold:
                lists['average'].append(data.get('average', None))
                lists['count'].append(data.get('count', None))
                lists['max'].append(data.get('max', None))
                lists['min'].append(data.get('min', None))
                lists['total'].append(data.get('total', None))
                lists['variance'].append(data.get('variance', None))
                lists['rel_frac'].append(round(relative_fractional_width * 100, 2))
                lists['abs_frac'].append(round(fractional_width * 100, 2))
                lists['level'].append(level)
                lists['timer'].append(key)
                lists['full-timer'].append(current_key_chain)
                fill_table_recursively(data['childs'], lists,
                                       start_location=current_location, level=level + 1,
                                       parent_fractional_width=fractional_width, parent_keys=current_key_chain,
                                       rel_threshold=rel_threshold, abs_threshold=abs_threshold, max_depth=max_depth)
            current_location += fractional_width
    return lists


def find_target_recursively(json_dict: dict, targets: str):
    for key in json_dict:
        data = json_dict[key]
        if isinstance(data, dict):
            if key in targets:
                return [True, data['childs']]
            [found, data] = find_target_recursively(data['childs'], targets)
            if found:
                return [True, data]
    return [False, json_dict]


def main():
    parser = argparse.ArgumentParser(description='Tool to plot walberla TimingTree JSON output files as pie chart.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('file', type=str,
                        help='JSON file generated from a walberla TimingTree')
    parser.add_argument('-r', '--rel-threshold', type=float, default=0.2,
                        help='threshold for the relative amount of time spent to be included in the plot')
    parser.add_argument('-a', '--abs-threshold', type=float, default=0.05,
                        help='threshold for the absolute amount of time spent to be included in the plot')
    outputs = parser.add_mutually_exclusive_group(required=False)
    outputs.add_argument('-o', '--output-file', type=str,
                         help='writes the pie chart to file')
    outputs.add_argument('--output-folder', type=str,
                         help='writes the pie chart to file')
    outputs.add_argument('--auto', action='store_true',
                        help='writes the pie chart to file')
    parser.add_argument('--max-depth', type=int, default=float('inf'),
                        help='maximal timing-tree depth to traverse')
    parser.add_argument('--target', type=str, default='', nargs='+',
                        help='starts traversing the timing tree at timer "target"')
    parser.add_argument('--show-filename', action='store_true',
                        help='displays the passed file name in the image')
    parser.add_argument('--preset', type=str, choices=['inner', 'outer', 'detail', 'vcycle', 'coarse', 'solve', 'kernel', 'smoother', 'fineLvl', 'cheby', 'apply', 'test'], default=None,
                        help='chooses one of three presets')
    parser.add_argument('--kernelLvl', type=int, default=6,
                        help='only relevant if preset=kernel, then it represents the matvec kernel on level <kernelLvl>')

    args = parser.parse_args()
    targets = ' '.join(args.target)

    if not os.path.exists(args.file):
        sys.exit('Input file does not exist.')

    lists = {
        'average': [],
        'count': [],
        'max': [],
        'min': [],
        'total': [],
        'variance': [],
        'rel_frac': [],
        'abs_frac': [],
        'level': [],
        'timer': [],
        'full-timer': [],
    }

    with open(args.file) as f:
        print(args.file)
        data = json.load(f)

        _, data = find_target_recursively(data, targets)
        lists = fill_table_recursively(data, lists,
                                       rel_threshold=args.rel_threshold, abs_threshold=args.abs_threshold,
                                       max_depth=args.max_depth)

    df = pd.DataFrame(lists)

    if args.output_file:
        with open(args.output_file + '.tex', 'w') as f:
            f.write(df.to_latex())
    elif args.auto or args.output_folder:
        input_path = Path(args.file)
        filename = input_path.stem
        out_path = str(input_path.parent) + "/../timings"
        if args.output_folder:
            out_path = args.output_folder
        Path(out_path).mkdir(parents=True, exist_ok=True)
        with open(out_path + '/' + filename + '.tex', 'w') as f:
            f.write(df.to_latex())
    else:
        # Display the DataFrame
        print(df)

    # Define the timers you want to filter
    desired_timers = [
        'Initialize Problem',                   ##
        'Initialization Outer solver',          ###
        'OuterLoop',                            ###
        'Iterative Refinement Solver',          ####      #* this is the main metric I want to reduce.
        'Inner Residual Computation',           #####
        'Copying-Casting',                      #####
        'InnerLoop',                            #####
        'Geometric Multigrid Solver',           ####      #* By finding the best trait-off between this and #V-cycles
        'Prolongation',                         #####
        'Restriction',                          #####
        'Smoother',                             #####
        'OuterResidual',                        ####
        'OuterError',                           ####
        'Compute error',                        ###
    ]
    if args.preset == 'outer':
        desired_timers = [
            'OuterLoop',                            ###
        ]
    elif args.preset == 'inner':
        desired_timers = [
            'Iterative Refinement Solver',          ####      #* this is the main metric I want to reduce.
            'Geometric Multigrid Solver',           ####      #* By finding the best trait-off between this and #V-cycles
            'OuterResidual',                        ####
            'OuterError',                           ####
        ]
    elif args.preset == 'detail':
        desired_timers = [
            'Inner Residual Computation',           #####
            'Copying-Casting',                      #####
            'InnerLoop',                            #####
        ]
    elif args.preset == 'vcycle':
        desired_timers = [
            'Geometric Multigrid Solver',           ####      #* By finding the best trait-off between this and #V-cycles
        ]
    elif args.preset == 'coarse':
        desired_timers = [
            'CG Solver',
        ]
    elif args.preset == 'solve':
        desired_timers = [
            'OuterSolve',
        ]
    elif args.preset == 'kernel':
        lvl = str(args.kernelLvl)
        desired_timers = [
            f'kernel{lvl}',
        ]
        df = df[df['full-timer'].str.contains('Smoother')]
    elif args.preset == 'smoother':
        lvl = str(args.kernelLvl)
        desired_timers = [
            'Smoother',
        ]
        df = df[df['full-timer'].str.endswith(f'Level {lvl}.Smoother')]
    elif args.preset == 'fineLvl':
        lvl = str(args.kernelLvl)
        desired_timers = [
            f'Level {lvl}',
        ]
        df = df[df['full-timer'].str.endswith(f'Geometric-Multigrid-Solver.Level-{lvl}')]
    elif args.preset == 'cheby':
        lvl = str(args.kernelLvl)
        desired_timers = [
            'apply',
            'assign',
            'elementwise multiply',
        ]
        suffixes = tuple(f'Level-{lvl}.Smoother.Chebyshev Smoother.{timer}' for timer in desired_timers)
        df = df[df['full-timer'].str.endswith(suffixes)]
    elif args.preset == 'apply':
        lvl = str(args.kernelLvl)
        desired_timers = [
            f'kernel{lvl}',
            'pre-communication',
            'post-communication',
        ]
        suffixes = [f'Level {lvl}.Smoother.Chebyshev Smoother.apply.Operator P1Function / VertexDoFFunction to P1Function / VertexDoFFunction.apply.{timer}' for timer in desired_timers]
        suffixes = tuple(suffixes + [f'Level {lvl}.Smoother.Operator P1Function / VertexDoFFunction to P1Function / VertexDoFFunction.apply.{timer}' for timer in desired_timers])
        df = df[df['full-timer'].str.endswith(suffixes)]
    elif args.preset == 'test':
        lvl = str(args.kernelLvl)
        # suffixes = tuple([f'Level {lvl}.Smoother.Chebyshev Smoother.apply.Operator P1Function / VertexDoFFunction to P1Function / VertexDoFFunction.apply.', f'Level {lvl}.Smoother.Operator P1Function / VertexDoFFunction to P1Function / VertexDoFFunction.apply.'])
        suffixes = f'Level {lvl}.Smoother.Chebyshev Smoother.apply.Operator P1Function / VertexDoFFunction to P1Function / VertexDoFFunction.apply.'
        df = df[df['full-timer'].str.contains(suffixes)]
    desired_timers = df['timer']

    # df = df[df['timer'].isin(desired_timers)]

    # Create a DataFrame with default values for the missing names
    default_values = {
        'timer': [timer for timer in desired_timers if timer not in df['timer'].unique()],
        'average': 0.,
        'count': -1,
        'max': 0.,
        'min': 0.,
        'total': 0.,
        'variance': 0.,
        'rel_frac': 0.,
        'abs_frac': 0.,
        'level': -1,
        'full-timer': 'This timer is not part of this algorithm',
    }
    # Concatenate the filtered DataFrame and the DataFrame with default values
    df = pd.concat([df, pd.DataFrame(default_values)])

    df = df.sort_values(by=['timer'])
    # Formatting
    # df['timer'] = df['timer'].apply(replace_whitespaces)
    # df['full-timer'] = df['full-timer'].apply(replace_whitespaces)
    # df['timer'] = df['timer'].astype(str).str.pad(max(df['timer'].astype(str).apply(len)), side='right', fillchar=' ')
    # df['full-timer'] = df['full-timer'].astype(str).str.pad(max(df['full-timer'].astype(str).apply(len)), side='right', fillchar=' ')
    # df['count'] = df['count'].astype(str).str.pad(max(df['count'].astype(str).apply(len)), side='left', fillchar=' ')
    # df['level'] = df['level'].astype(str).str.pad(max(df['level'].astype(str).apply(len)), side='left', fillchar=' ')
    df.style \
        .format(precision=3, thousands="'", decimal=".") \
        .format_index(str.upper, axis=1)
    df.set_index('timer', inplace=True)

    if args.output_file:
        df.to_csv(args.output_file + '.dat', sep='\t', na_rep='', float_format='%.4E', header=True)
    elif args.auto or args.output_folder:
        input_path = Path(args.file)
        filename = input_path.stem
        out_path = str(input_path.parent) + "/../timings"
        if args.output_folder:
            out_path = args.output_folder
        Path(out_path).mkdir(parents=True, exist_ok=True)
        df.to_csv(out_path + '/' + filename + '.dat', sep='\t', na_rep='', float_format='%.4E', header=True)
    else:
        # Display the DataFrame
        print(df)


if __name__ == '__main__':
    main()
