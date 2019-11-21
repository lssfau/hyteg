import csv
import collections
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import argparse


def main(datafile, perfgroup):
    for messure in [perfgroup]:
        likwidRegions = collections.OrderedDict()

        with open(datafile, 'r') as csvfile:
            lines = csv.reader(csvfile, delimiter=',')
            region = 'walberla'
            likwidRegions[region] = {}
            for l in lines:
                if (l[0] == 'Region' or l[0] == 'TABLE'):
                    region = l[1]
                    likwidRegions[region] = {}
                else:
                    likwidRegions[region][l[0]] = l[1:]

        del likwidRegions['walberla']

        numberOfProcessors = []
        outputData = collections.defaultdict(list)

        print(likwidRegions.keys())
        for region in likwidRegions:
            if 'Region calls STAT' not in likwidRegions[region]:
                procs = 1
            else:
                procs = int(likwidRegions[region]['Region calls STAT'][1]) / \
                    int(likwidRegions[region]['Region calls STAT'][2])
            numberOfProcessors.append(procs)

            outputData[region.replace("Region ", "").split("-level")[0]].append(
                float(likwidRegions[region][messure][0]))

        totalProcs = list(set(numberOfProcessors))

        for region in outputData:
            plt.plot(totalProcs, outputData[region], '--o', label=region)
        plt.xlabel('Procs')

        plt.ylabel(messure)

        plt.ylim(bottom=0)
        plt.xlim(left=0, right=totalProcs[-1] + 1)
        plt.xticks(range(int(totalProcs[0]+1), int(totalProcs[-1] + 1), 2))
        plt.grid(True)
        plt.gca().get_yaxis().set_major_formatter(
            ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

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
    # parser.add_argument("minLevel", help="first level to plot", type=int)
    # parser.add_argument("maxLevel", help="last level to plot", type=int)
    parser.add_argument("--perfgroup", default='MFLOP/s STAT',
                        help="which likwid performance group to plot")
    args = parser.parse_args()

    main(args.datafile, args.perfgroup)
