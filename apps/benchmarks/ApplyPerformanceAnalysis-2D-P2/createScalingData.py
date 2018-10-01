import csv
import subprocess
import socket
import collections
#import matplotlib.pyplot as plt


def main():

    #procs = range(20)
    totalProcs = range(2,20)

    if "emmy" in socket.gethostname():
        subprocess.call(["source", "~/script/modules/intel.emmy.modules"])

    f = open('createScaling.txt','w')

    for i in totalProcs:
        subprocess.call(["likwid-mpirun", "-np", str(i), "-g", "MEM", "-m", "-O",
                         "./ApplyPerformanceAnalysis-2D-P2", "ApplyPerformanceAnalysis-2D-P2.prm",
                         "-Parameters.level=11"], stdout=f)
    f.close()

if __name__ == "__main__":
    main()
