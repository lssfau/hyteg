import csv
import subprocess
import socket
import collections
#import matplotlib.pyplot as plt


def main():
    totalProcs = range(2,41)

    f = open('scalingData.txt','w')
    for i in totalProcs:
        #if i is 33:
            continue
        subprocess.run(["likwid-mpirun", "-np", str(i), "-g", "MEM_DP", "-m", "-O","./ApplyPerformanceAnalysis-2D-P2", "ApplyPerformanceAnalysis-2D-P2.prm", "-Parameters.level=11"], stdout=f)
    f.close()

if __name__ == "__main__":
    main()
