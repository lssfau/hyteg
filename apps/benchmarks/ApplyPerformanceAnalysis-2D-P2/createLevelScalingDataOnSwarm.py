import subprocess
import os


def main(cores, minLevel, maxLevel):
    filename = os.environ["DOCKER_IMAGE_NAME"]
    f = open(filename + '.txt', 'w')

    for level in range(minLevel, maxLevel + 1):
        levelstring = str("-Parameters.level=" + str(level))
        subprocess.run(["likwid-perfctr", "-C", "1", "-g", "MEM_DP", "-m", "-O", "--stats", "./ApplyPerformanceAnalysis-2D-P2",
                            "ApplyPerformanceAnalysis-2D-P2.prm", levelstring], stdout=f)

    f.close()


if __name__ == "__main__":
    main(1, 2, 15)
