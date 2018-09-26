import os
import time
import math
import random
import csv
import re
from influxdb import InfluxDBClient
from git import Repo


def main():
    # try:
    #     write_user_pw = os.environ["INFLUXDB_WRITE_USER"]
    # except KeyError:
    #     import sys
    #     print('Password for the InfluxDB write_user was not set.\n',
    #           'See https://docs.gitlab.com/ee/ci/variables/#secret-variables', file=sys.stderr)
    #     exc_info = sys.exc_info()
    #     raise exc_info[0].with_traceback(exc_info[1], exc_info[2])
    write_user_pw = 'hyteg'

    client = InfluxDBClient('i10grafana.informatik.uni-erlangen.de', 8086,
                            'terraneo', write_user_pw, 'terraneo')

    #repo = Repo(search_parent_directories=True)
    #commit = repo.head.commit

    likwidRegions = {}

    with open('./output.txt', 'r') as csvfile:
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

    for region in likwidRegions:
        procs = int(likwidRegions[region]['Region calls STAT'][1]) / int(likwidRegions[region]['Region calls STAT'][2])

        json_body += [
            {
                'measurement': 'P2_Apply_Benchmark',
                'tags': {
                    'host': os.uname()[1],
                    'project': 'terraneo',
                    #'image': os.environ["DOCKER_IMAGE_NAME"],
                    'Processes': int(procs),
                    #'commit': commit,
                    'region': region,
                },
                'time': int(time.time()) - int(procs),
                'fields': {'DP MFLOP/s': float(likwidRegions[region]['DP MFLOP/s STAT'][0])}
                           #'Memory bandwidth [MBytes/s]': float(likwidRegions[region]['Memory bandwidth [MBytes/s] STAT'][0])}
            }
        ]

    print(json_body)
    client.write_points(json_body, time_precision='s')


if __name__ == "__main__":
    main()
