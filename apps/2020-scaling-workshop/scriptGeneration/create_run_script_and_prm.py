import argparse
import time
from templates import job_template, param_template


def get_script(variables):
    return job_template.format(**variables)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("num_nodes", help="number of nodes to be used", type=int)
    parser.add_argument("numEPS", help="number of edges per side", type=int)
    parser.add_argument("--maxLevel", default=8,
                        help="max level for multigird. ATTENTION: P2 is reduced by 1", type=int)
    parser.add_argument("--num_cores", default=32,
                        help="number of cores per node", type=int)
    parser.add_argument("--walltime", default="00:10:00",
                        help="walltime for the job")

    args = parser.parse_args()
    args_dict = vars(args)
    args_dict['total_num_procs'] = args_dict['num_nodes'] * args_dict['num_cores']
    datestamp = time.strftime("%y_%m_%d-%H.%M.%S")
    total_tets = args_dict['numEPS'] ** 3 * 24
    print("{} cores; {} nodes; {} tets".format(args_dict['num_cores'], args_dict['num_nodes'], total_tets))
    print("{} tets per node\n{} tets per core".format(total_tets / args_dict['num_nodes'], total_tets / args_dict['total_num_procs']))
    args_dict['job_name'] = "{}-FMG-Scaling-{}nodes-{}cores-maxLevel{}-{}tets".format(datestamp, args_dict['num_nodes'],
                                                                                      args_dict['maxLevel'],
                                                                                      args_dict['num_cores'], total_tets)
    args_dict['paramfile_name'] = args_dict['job_name'] + '.prm'

    with open("../hawk/" + args_dict['job_name'] + '.job', 'w') as f:
        f.write(job_template.format(**args_dict))

    with open("../hawk/" + args_dict['job_name'] + '.prm', 'w') as f:
        f.write(param_template.format(**args_dict))
    # print(param_template.format(**args_dict))
    # print(job_template.format(**args_dict))
