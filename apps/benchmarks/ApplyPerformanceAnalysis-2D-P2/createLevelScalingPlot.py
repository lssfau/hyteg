import csv
import collections
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import argparse


def main(datafile,outputfile,minLevel,maxLevel):
    for messure in ['MFLOP/s STAT']:
        likwidRegions = collections.OrderedDict()

        with open(datafile, 'r') as csvfile:
            lines = csv.reader(csvfile, delimiter=',')
            region = 'walberla'
            likwidRegions[region] = {}
            for l in lines:
                if (l[0] == 'TABLE'):
                    region = l[1]
                    likwidRegions[region] = {}
                else:
                    likwidRegions[region][l[0]] = l[1:]

        del likwidRegions['walberla']

        levelRange = range(minLevel,maxLevel+1)
        outputData = collections.defaultdict(list)

        for region in likwidRegions:
            outputData[region.split("-level")[0]].append(float(likwidRegions[region][messure][0]))

        print(outputData)

        for region in outputData:
            plt.plot(levelRange, outputData[region][minLevel-2:], '--o', label=region)
        plt.xlabel('Refinement Level')

        plt.ylabel("Performance (Mflops/s)")

        plt.ylim(bottom=0)
        # plt.xlim(left=0, right=totalProcs[-1] + 1)
        # plt.xticks(range(2, int(totalProcs[-1] + 1), 2))
        plt.grid(True)
        plt.gca().get_yaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

        plt.legend()
        namepart = messure.split('[')[0].replace(' ', '_').replace('/', '-')
        savename = namepart + outputfile
        plt.tight_layout()
        plt.savefig(savename)
        plt.clf()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("datafile", help="file to read input from")
    parser.add_argument("outputfile", help="filename for output including the ending; something will be added")
    parser.add_argument("minLevel", help="first level to plot", type=int)
    parser.add_argument("maxLevel", help="last level to plot", type=int)
    args = parser.parse_args()

    main(args.datafile,args.outputfile,args.minLevel,args.maxLevel)
