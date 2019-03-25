import subprocess
import time
import argparse


def main(cores, minLevel, maxLevel):
    datestamp = time.strftime("%y_%m_%d-%I:%M:%S")
    f = open(datestamp + 'levelScalingData.txt', 'w')

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
    parser.add_argument("cores", help="number of cores to run with")
    parser.add_argument("minLevel", default=2, help="minimal scaling Level")
    parser.add_argument("maxLevel", default=15, help="maximal scaling Level")
    args = parser.parse_args()
    main()
