import os
import time
import math
import random
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

    repo = Repo(search_parent_directories=True)
    commit = repo.head.commit

    with open("P1BenchmarkOutput.txt") as f:
        s = f.read()
    apply = re.search('apply runtime: (\d*.\d*)', s)
    assign = re.search('assign runtime: (\d*.\d*)', s)
    interpolate = re.search('interpolate runtime: (\d*.\d*)', s)
    l = re.search('level: (\d*.\d*)', s)

    json_body = [
        {
            'measurement': 'P1_Benchmark',
            'tags': {
                'host'     : os.uname()[1],
                'project'  : 'terraneo',
                'image'    : os.environ["DOCKER_IMAGE_NAME"],
                'Level'    : int(l.group(1)),
                'commit'  : commit,
            },
            'time': int(time.time()),
            'fields': {'apply_runtime': float(apply.group(1)),
                       'assign_runtime': float(assign.group(1)),
                       'interpolate_runtime': float(interpolate.group(1))}
        }
    ]
    client.write_points(json_body, time_precision='s')


if __name__ == "__main__":
    main()
