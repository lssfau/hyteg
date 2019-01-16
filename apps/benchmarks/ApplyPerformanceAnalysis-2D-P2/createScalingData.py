import csv
import subprocess
import socket
import collections
#import matplotlib.pyplot as plt


def main():

    #procs = range(20)
    totalProcs = range(2,20)

    # if "emmy" in socket.gethostname():
    #     subprocess.call(["source", "~/script/modules/gcc.emmy.modules"])

    f = open('scalingData.txt','w')

    for i in totalProcs:
        subprocess.call(["likwid-mpirun", "-np", str(i), "-g", "MEM_DP", "-m", "-O",
                         "./ApplyPerformanceAnalysis-2D-P2", "ApplyPerformanceAnalysis-2D-P2.prm",
                         "-Parameters.level=11"], stdout=f)
    f.close()

if __name__ == "__main__":
    main()
