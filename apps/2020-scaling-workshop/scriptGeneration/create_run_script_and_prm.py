import argparse
import time
import sys
import os

from templates import parameter_file_01_cube, job_file_hawk, job_file_supermuc


def create_files(args, args_dict):
    datestamp = time.strftime("%y_%m_%d-%H.%M.%S")
    args_dict['prefix'] = datestamp

    # parameter file
    if args.benchmark == 'cube':
        binary_name = 'Scaling_Benchmark_01_Cube'
        num_tets = (args.num_edges_per_side ** 3) * 24

        base_name = '_'.join([datestamp, 'benchmark_cube', f'{args.num_nodes}_nodes', f'{num_tets}_tets'])
        parameter_file_name = base_name + '.prm'
        args_dict['db_file'] = base_name + '.db'

        parameter_file = parameter_file_01_cube(**args_dict)
    else:
        print('Invalid benchmark. Specify subargument!')
        sys.exit(1)

    # job file
    if args.machine == 'hawk':
        job_file_name = base_name + '.job'
        args_dict['out_dir'] = '../hawk'
        args_dict['binary_name'] = binary_name
        args_dict['job_name'] = base_name
        args_dict['paramfile_name'] = parameter_file_name
        args_dict['total_num_procs'] = args.num_cores * args.num_nodes
        job_file = job_file_hawk(**args_dict)
    elif args.machine == 'supermuc':
        job_file_name = base_name + '.job'
        args_dict['out_dir'] = '../supermuc'
        args_dict['binary_name'] = binary_name
        args_dict['job_name'] = base_name
        args_dict['paramfile_name'] = parameter_file_name
        job_file = job_file_supermuc(**args_dict)
    else:
        print('Invalid machine.')
        sys.exit(1)

    with open(os.path.join(args.out_dir, parameter_file_name), 'w') as f:
        f.write(parameter_file)

    with open(os.path.join(args.out_dir, job_file_name), 'w') as f:
        f.write(job_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--machine", help="hawk or supermuc", type=str, required=True)
    parser.add_argument("--num_nodes", help="number of nodes to be used", type=int, required=True)
    parser.add_argument("--num_cores", help="number of cores per node", type=int, required=True)
    parser.add_argument("--walltime", default="00:10:00", help="walltime for the job")

    subparsers = parser.add_subparsers(dest='benchmark')

    parser_cube = subparsers.add_parser('cube', help='cube benchmark')

    parser_cube.add_argument("--scenario", type=int, required=True)
    parser_cube.add_argument("--num_edges_per_side", help="number of edges per side", type=int, required=True)
    parser_cube.add_argument("--max_level", default=7, help="max level for multigrid", type=int, required=True)

    args = parser.parse_args()
    args_dict = vars(args)

    create_files(args, args_dict)
