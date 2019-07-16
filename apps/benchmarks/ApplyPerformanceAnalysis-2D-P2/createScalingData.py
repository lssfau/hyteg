import subprocess
import argparse
import time


def processorScaling(minCores, maxCores, level, outputfile):
    datestamp = time.strftime("%y_%m_%d-%I.%M.%S")
    f = open(datestamp + outputfile, 'w')

    for i in range(minCores, maxCores + 1):
        if i is 1:
            subprocess.run(
                ["likwid-mpirun", "-np", "1", "likwid-perfctr", "-g", "FLOPS_DP", "-m", "-O", "--stats", "./ApplyPerformanceAnalysis-2D-P2",
                 "ApplyPerformanceAnalysis-2D-P2.prm", f"-Parameters.level={level}"], stdout=f)
        else:
            subprocess.run(["likwid-mpirun", "-np", f"{i}", "-g", "FLOPS_DP", "-m", "-O", "./ApplyPerformanceAnalysis-2D-P2",
                            "ApplyPerformanceAnalysis-2D-P2.prm", f"-Parameters.level={level}"], stdout=f)

    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("minCores", help="minimal number of cores", type=int)
    parser.add_argument("maxCores", help="maximal number of cores", type=int)
    parser.add_argument("level", help="refinement level", type=int)
    parser.add_argument("outputfile", help="filename for output including the ending; something will be added")
    args = parser.parse_args()
    processorScaling(args.minCores, args.maxCores, args.level, args.outputfile)
