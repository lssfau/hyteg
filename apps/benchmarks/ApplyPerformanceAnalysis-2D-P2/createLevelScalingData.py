import subprocess
import time
import argparse


def main(cores, minLevel, maxLevel, outputfile):
    datestamp = time.strftime("%y_%m_%d-%I:%M:%S")
    f = open(datestamp + outputfile, 'w')

    if cores is 1:
        for level in range(minLevel, maxLevel + 1):
            levelstring =  str("-Parameters.level=" + str(level))
            subprocess.run(["likwid-perfctr", "-C", "1", "-g", "MEM_DP", "-m", "-O", "--stats", "./ApplyPerformanceAnalysis-2D-P2",
            "ApplyPerformanceAnalysis-2D-P2.prm", levelstring], stdout=f)
    else:
        for level in range(minLevel, maxLevel + 1):
            levelstring =  str("-Parameters.level=" + str(level))
            subprocess.run(["likwid-mpirun", "-n", str(cores), "-g", "MEM_DP", "-m", "-O", "./ApplyPerformanceAnalysis-2D-P2",
                            "ApplyPerformanceAnalysis-2D-P2.prm", levelstring], stdout=f)
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cores", help="number of cores to run with", type=int)
    parser.add_argument("minLevel", help="minimal scaling Level", type=int)
    parser.add_argument("maxLevel", help="maximal scaling Level", type=int)
    parser.add_argument("outputfile", help="filename for output including the ending; something will be added")
    args = parser.parse_args()
    main(args.cores, args.minLevel, args.maxLevel, args.outputfile)
