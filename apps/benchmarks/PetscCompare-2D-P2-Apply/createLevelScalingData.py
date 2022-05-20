import subprocess
import time
import argparse


def main(cores, minLevel, maxLevel, outputfile):
    datestamp = time.strftime("%y_%m_%d-%I.%M.%S")
    f = open(datestamp + outputfile, 'w')

    if cores == 1:
        for level in range(minLevel, maxLevel + 1):
            levelstring = str("-Parameters.level=" + str(level))
            subprocess.run(["likwid-mpirun", "-n", "1", "likwid-perfctr", "-g", "MEM_DP", "-m", "-O", "--stats", "./PetscCompare-2D-P2-Apply.cpp",
                            "PetscCompare-2D-P2-Apply.prm", levelstring], stdout=f)
    else:
        for level in range(minLevel, maxLevel + 1):
            levelstring = str("-Parameters.level=" + str(level))
            subprocess.run(["likwid-mpirun", "-n", str(cores), "-g", "MEM_DP", "-m", "-O", "./PetscCompare-2D-P2-Apply",
                            "PetscCompare-2D-P2-Apply.prm", levelstring], stdout=f)
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cores", help="number of cores to run with", type=int)
    parser.add_argument(
        "outputfile", help="filename for output including the ending; something will be added")
    parser.add_argument("--minLevel", default=2,
                        help="minimal scaling Level", type=int)
    parser.add_argument("--maxLevel", default=15,
                        help="maximal scaling Level", type=int)
    args = parser.parse_args()
    main(args.cores, args.minLevel, args.maxLevel, args.outputfile)
