import argparse
import time
import sys
import os

from templates import create_parameter_file, job_file_hawk, job_file_supermuc


def create_files(args, args_dict):
    datestamp = time.strftime('%y_%m_%d-%H_%M_%S')

    # parameter file
    binary_name = 'Benchmark_05_PipeScaling'

    base_name = '_'.join(['moc', datestamp, f'nodes_{args.num_nodes}', f'tets_{args.diameter_cubes**2 * args.length_cubes * 6}'])
    parameter_file_name = base_name + '.prm'
    args_dict['db_file_name'] = base_name + '.db'

    parameter_file = create_parameter_file(**args_dict)

    # job file
    if args.machine == 'hawk':
        raise Exception('Hawk untested')
        job_file_name = base_name + '.job'
        args_dict['out_dir'] = '../hawk'
        args_dict['path'] = os.getcwd().rstrip('scriptGeneration')
        args_dict['binary_name'] = binary_name
        args_dict['job_name'] = base_name
        args_dict['paramfile_name'] = parameter_file_name
        args_dict['total_num_procs'] = args.num_cores * args.num_nodes
        job_file = job_file_hawk(**args_dict)
    elif args.machine == 'supermuc':
        job_file_name = base_name + '.job'
        args_dict['out_dir'] = '.'
        args_dict['binary_name'] = binary_name
        args_dict['job_name'] = base_name
        args_dict['paramfile_name'] = parameter_file_name
        job_file = job_file_supermuc(**args_dict)
    else:
        print('Invalid machine.')
        sys.exit(1)

    parameter_file_output_path = os.path.join(args.out_dir, parameter_file_name)
    with open(parameter_file_output_path, 'w') as f:
        f.write(parameter_file)
        print(f'Written .prm file to {parameter_file_output_path}')

    job_file_output_path = os.path.join(args.out_dir, job_file_name)
    with open(job_file_output_path, 'w') as f:
        f.write(job_file)
        print(f'Written .job file to {job_file_output_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--machine', help='hawk or supermuc', type=str, required=True)
    parser.add_argument('--num_nodes', help='number of nodes to be used', type=int, required=True)
    parser.add_argument('--num_mpi_procs_per_node', help='number of MPI processes per node', type=int, required=True)
    parser.add_argument('--walltime', default='00:10:00', help='walltime for the job')

    parser.add_argument('--level', help='refinement level', type=int, required=True)
    parser.add_argument('--diameter_cubes', help='number of cubes in y and z direction', type=int, required=True)
    parser.add_argument('--length_cubes', help='number of cubes in x direction', type=int, required=True)
    parser.add_argument('--num_time_steps', default=10, help='number of time steps', type=int )

    parser.add_argument('--lb_type', default=0, help='load balancing algorithm', type=int )
    parser.add_argument('--lb_part_x_size', default=1.0, help='load balancing x cube size', type=float )

    args = parser.parse_args()
    args_dict = vars(args)

    create_files(args, args_dict)
