import csv
import subprocess
import socket
import collections
import matplotlib.pyplot as plt


def main():
    totalProcs = range(2,20)
    likwidRegions = collections.OrderedDict()

    with open('scalingData.txt','r') as csvfile:
        lines = csv.reader(csvfile, delimiter=',')
        region = 'walberla'
        likwidRegions[region] = {}
        for l in lines:
            if(l[0] == 'Region'):
                region = l[1]
                likwidRegions[region] = {}
            else:
                likwidRegions[region][l[0]] = l[1:]

    json_body = []

    del likwidRegions['walberla']

    #print(likwidRegions)
    #client.write_points(json_body, time_precision='s')


    numberOfProcessors = []
    outputData = collections.defaultdict(list)

    for region in likwidRegions:
        print(region)
        procs = int(likwidRegions[region]['Region calls STAT'][1]) / int(likwidRegions[region]['Region calls STAT'][2])
        numberOfProcessors.append(procs)

        outputData[region.split("-level")[0]].append(float(likwidRegions[region]['Memory bandwidth [MBytes/s] STAT'][0]))

        # json_body += [
        #     {
        #         'measurement': 'P2_Apply_Benchmark',
        #         'tags': {
        #             'host': os.uname()[1],
        #             'project': 'terraneo',
        #             #'image': os.environ["DOCKER_IMAGE_NAME"],
        #             'Processes': int(procs),
        #             #'commit': commit,
        #             'region': region,
        #         },
        #         'time': int(time.time()) - int(procs),
        #         'fields': {'DP MFLOP/s': float(likwidRegions[region]['DP MFLOP/s STAT'][0])}
        #                    #'Memory bandwidth [MBytes/s]': float(likwidRegions[region]['Memory bandwidth [MBytes/s] STAT'][0])}
        #     }
        # ]

    print(numberOfProcessors)
    print(outputData)
    #client.write_points(json_body, time_precision='s')

    for region in outputData:
        print(region)
        print(totalProcs)
        print(outputData[region])
        plt.plot(totalProcs,outputData[region], '--o', label=region)
    # plt.ylabel('petsc / hyteg (speedup)')
    plt.xlabel('Procs')
    plt.ylabel('Bandwidth [MByte/s]')
    plt.ylim(bottom=0)
    plt.xlim(left=0,right=totalProcs[-1]+1)
    plt.xticks(totalProcs)

    plt.legend()
    savename = 'scaling.png'
    plt.savefig(savename)

    #plt.show()



if __name__ == "__main__":
    main()
