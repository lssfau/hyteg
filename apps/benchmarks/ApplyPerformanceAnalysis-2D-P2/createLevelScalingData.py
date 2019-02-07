import subprocess


def main():
    alllevel = range(2, 14)

    f = open('levelScalingData.txt', 'w')
    for level in alllevel:
        subprocess.run(["likwid-perfctr", "-C", 1, "-g", "MEM_DP", "-m", "-O", "./ApplyPerformanceAnalysis-2D-P2",
                        "ApplyPerformanceAnalysis-2D-P2.prm", "-Parameters.level=" + level], stdout=f)
    f.close()


if __name__ == "__main__":
    main()
