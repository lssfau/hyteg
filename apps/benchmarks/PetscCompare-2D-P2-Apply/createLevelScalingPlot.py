import csv
import collections
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import argparse


def main(datafile, perfgroup, minLevel, maxLevel):
    # ecmPredictionFixed = [12.04, 12.04, 12.04, 12.04, 12.04, 8.46, 7.54, 7.26, 6.44]
    # ecmPrediction = [12.04, 12.04, 12.04, 12.04, 12.04, 7.26, 7.26, 7.26, 6.44]
    # ecmPrediction = [i * 1000 for i in ecmPrediction]

    for messure in [perfgroup]:
        likwidRegions = collections.OrderedDict()

        with open(datafile, 'r') as csvfile:
            lines = csv.reader(csvfile, delimiter=',')
            region = 'walberla'
            likwidRegions[region] = {}
            for l in lines:
                if not l:
                    continue
                if l[0] == 'TABLE':
                    region = l[1]
                    likwidRegions[region] = {}
                else:
                    likwidRegions[region][l[0]] = l[1:]

        del likwidRegions['walberla']

        levelRange = range(minLevel, maxLevel + 1)
        print("level Range:")
        print(levelRange)
        outputData = collections.defaultdict(list)

        max_measure_level = 0
        min_measure_level = 100
        for region in likwidRegions:
            if "-level" not in region:
                continue
            current_level = int(region.split("-level")[1][0:2])
            print(current_level)
            if max_measure_level < current_level:
                max_measure_level = current_level
            if min_measure_level > current_level:
                min_measure_level = current_level
            outputData[region.split(
                "-level")[0]].append(float(likwidRegions[region][messure][0]))

        print("Max measure level: " + str(max_measure_level))
        print("Min measure level: " + str(min_measure_level))

        print(outputData)

        # plt.plot(levelRange, ecmPrediction, '--', label="ecm prediction v to v", color='black')

        # for region in outputData:
        plt.plot(levelRange, outputData["Region HyTeG-apply"][minLevel - min_measure_level:maxLevel + 1 - min_measure_level], '--o',
                 label="Region HyTeG-apply")
        plt.plot(levelRange, outputData["Region Petsc-MatMult"][minLevel - min_measure_level:maxLevel + 1 - min_measure_level], '--o',
                 label="Region Petsc-MatMult")

        plt.xlabel('Refinement Level')

        plt.ylabel("Performance (Mflops/s)")

        # plt.ylim(bottom=0)
        # plt.xlim(left=0, right=totalProcs[-1] + 1)
        plt.xticks(range(minLevel, maxLevel+1))
        plt.grid(True)
        # plt.gca().get_yaxis().set_major_formatter(
        #     ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

        plt.legend()
        namepart = messure.split('[')[0].replace(' ', '_').replace('/', '-')
        savename = ".".join(datafile.split(
            ".")[0:-1]) + "-" + namepart + ".pdf"
        plt.tight_layout()
        plt.savefig(savename)
        plt.clf()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "datafile", help="file to read input from and also name of output")
    parser.add_argument("minLevel", help="first level to plot", type=int)
    parser.add_argument("maxLevel", help="last level to plot", type=int)
    parser.add_argument("--perfgroup", default='MFLOP/s STAT',
                        help="which likwid performance group to plot")
    args = parser.parse_args()

    main(args.datafile, args.perfgroup, args.minLevel, args.maxLevel)
