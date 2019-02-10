import csv
import collections
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def main():
    for messure in ['MFLOP/s STAT']:
        likwidRegions = collections.OrderedDict()

        with open('levelScalingData.txt', 'r') as csvfile:
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

        levelRange = range(2,15)
        outputData = collections.defaultdict(list)

        for region in likwidRegions:
            outputData[region.split("-level")[0]].append(float(likwidRegions[region][messure][0]))

        print(outputData)

        for region in outputData:
            plt.plot(levelRange, outputData[region], '--o', label=region)
        plt.xlabel('Procs')

        plt.ylabel(messure)

        plt.ylim(bottom=0)
        # plt.xlim(left=0, right=totalProcs[-1] + 1)
        # plt.xticks(range(2, int(totalProcs[-1] + 1), 2))
        plt.grid(True)
        plt.gca().get_yaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

        plt.legend()
        namepart = messure.split('[')[0].replace(' ', '_').replace('/', '-')
        savename = 'scaling_' + namepart + '.pdf'
        plt.tight_layout()
        plt.savefig(savename)
        plt.clf()


if __name__ == "__main__":
    main()
