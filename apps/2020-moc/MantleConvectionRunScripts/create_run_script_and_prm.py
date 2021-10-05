import argparse
import time
import sys
import os

from templates import create_parameter_file, job_file_hawk, job_file_supermuc


def create_files(args, args_dict):
    datestamp = time.strftime('%y_%m_%d-%H_%M_%S')

    # parameter file
    binary_name = 'MantleConvection'

    base_name = '_'.join(['mc', datestamp, f'nodes_{args.num_nodes}'])
    parameter_file_name = base_name + '.prm'
    args_dict['base_name'] = base_name
    args_dict['output_directory'] = args.out_dir

    parameter_file = create_parameter_file(**args_dict)

    # job file
    if args.machine == 'hawk':
        job_file_name = base_name + '.job'
        args_dict['out_dir'] = '.'
        args_dict['binary_name'] = binary_name
        args_dict['job_name'] = base_name
        args_dict['paramfile_name'] = parameter_file_name
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
    parser.add_argument('--walltime', default='01:00:00', help='walltime for the job')

    parser.add_argument('--ntan', default=5, help='parameter for spherical shell mesh', type=int)
    parser.add_argument('--nrad', default=4, help='parameter for spherical shell mesh', type=int)

    parser.add_argument('--max_level', default=4, help='max level for multigrid', type=int)
    parser.add_argument('--ra', default=1e8, help='Rayleigh number')
    parser.add_argument('--cfl', default=1.0, help='CFL number')
    parser.add_argument('--out_dir', default='.', help='output directory')
    parser.add_argument('--max_num_time_steps', default=10000, type=int)
    parser.add_argument('--uzawa_omega', default=0.3, type=float)

    parser.add_argument('--uzawa_pre', default=10, type=int)
    parser.add_argument('--uzawa_post', default=10, type=int)
    parser.add_argument('--uzawa_inner', default=6, type=int)

    parser.add_argument('--stokes_abs_tol', default=1e-8, type=float)
    parser.add_argument('--stokes_rel_tol', default=1e-8, type=float)

    parser.add_argument('--vtk_interval', default=10, type=int)
    parser.add_argument('--vtk_vertex_dofs', default=False, type=bool)

    parser.add_argument('--sph_tmp_interval', default=10, type=int)

    args = parser.parse_args()
    args_dict = vars(args)

    create_files(args, args_dict)
