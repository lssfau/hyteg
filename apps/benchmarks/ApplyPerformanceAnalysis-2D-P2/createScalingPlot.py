import csv
import collections
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def main():
    for messure in ['Memory bandwidth [MBytes/s] STAT','MFLOP/s STAT']:
        likwidRegions = collections.OrderedDict()

        with open('scalingData.txt', 'r') as csvfile:
            lines = csv.reader(csvfile, delimiter=',')
            region = 'walberla'
            likwidRegions[region] = {}
            for l in lines:
                if (l[0] == 'Region'):
                    region = l[1]
                    likwidRegions[region] = {}
                else:
                    likwidRegions[region][l[0]] = l[1:]

        del likwidRegions['walberla']

        numberOfProcessors = []
        outputData = collections.defaultdict(list)

        for region in likwidRegions:
            procs = int(likwidRegions[region]['Region calls STAT'][1]) / int(likwidRegions[region]['Region calls STAT'][2])
            numberOfProcessors.append(procs)

            outputData[region.split("-level")[0]].append(
                float(likwidRegions[region][messure][0]))

        totalProcs = list(set(numberOfProcessors))

        for region in outputData:
            plt.plot(totalProcs, outputData[region], '--o', label=region)
        plt.xlabel('Procs')

        plt.ylabel(messure)

        plt.ylim(bottom=0)
        plt.xlim(left=0, right=totalProcs[-1] + 1)
        plt.xticks(range(2, int(totalProcs[-1] + 1), 2))
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
