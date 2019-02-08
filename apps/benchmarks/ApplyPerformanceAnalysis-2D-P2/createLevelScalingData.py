import subprocess


def main():
    alllevel = range(2, 15)

    f = open('levelScalingData.txt', 'w')
    for level in alllevel:
        levelstring =  str("-Parameters.level=" + str(level))
        subprocess.run(["likwid-perfctr", "-C", "1", "-g", "MEM_DP", "-m", "-O", "--stats", "./ApplyPerformanceAnalysis-2D-P2",
                        "ApplyPerformanceAnalysis-2D-P2.prm", levelstring], stdout=f)
    f.close()


if __name__ == "__main__":
    main()
