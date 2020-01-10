import argparse
import subprocess

parser = argparse.ArgumentParser(description='Discover if the error scaling is correct.')
parser.add_argument('--max-level', help='The finest level we use.', type=int, default=6)
parser.add_argument('--num-processes', help='The number of processes we use for a calculation.', type=int, default=8)
parser.add_argument('--smoothing-steps', help='number of smoothing steps', type=int, default=2)
parser.add_argument('--coarse-iter', help='The number of coarse grid iterations.', type=int, default=15)

pargs = parser.parse_args()

smoother_order_pairs = [
    ("richardson", 0),
    ("jacobi", 0),
    ("chebyshev", 1),
    ("chebyshev", 2),
    ("chebyshev", 3),
    ("chebyshev", 4),
    ("chebyshev", 5),
]

for (smootherType, order) in smoother_order_pairs:
    print('!', smootherType, order)
    exe_name = 'P1GMGSmootherComparison'
    args = [
        'mpirun', '-np', '{}'.format(pargs.num_processes),
        './{}'.format(exe_name),
        '{}.prm'.format(exe_name),
        '-Parameters.maxLevel={}'.format(pargs.max_level),
        '-Parameters.writeCSV=true',
        '-Parameters.smootherType={}'.format(smootherType),
        '-Parameters.order={}'.format(order),
        '-Parameters.maxCoarseIter={}'.format(pargs.coarse_iter),
        '-Parameters.smoothingSteps={}'.format(pargs.smoothing_steps)]
    p = subprocess.Popen(args)
    p.communicate()
