import os
import time
import math
import random
import re
from influxdb import InfluxDBClient


def main():
    write_user_pw = 'hyteg'

    client = InfluxDBClient('i10grafana.informatik.uni-erlangen.de', 8086,
                            'terraneo', write_user_pw, 'terraneo')

    commit = os.environ["CI_COMMIT_SHA"]
    commit_short = os.environ["CI_COMMIT_SHORT_SHA"]
    branch = os.environ["CI_COMMIT_REF_NAME"]

    with open("./BuildTiming.txt") as f:
        s = f.read()
    hyteg = re.search('hyteg buildtime (\d*.\d*)', s)
    tests = re.search('tests buildtime (\d*.\d*)', s)
    apps = re.search('apps buildtime (\d*.\d*)', s)

    print(hyteg)

    json_body = [
        {
            'measurement': 'Build Timing',
            'tags': {
                'host': os.uname()[1],
                'project': 'terraneo',
                'commit': commit,
                'commit_short': commit_short,
                'branch': branch,
            },
            'time': int(time.time()),
            'fields': {'hyteg_buildtime': float(hyteg.group(1)),
                       'tests_buildtime': float(tests.group(1)),
                       'apps_buildtime': float(apps.group(1))}
        }
    ]
    print(json_body)
    client.write_points(json_body, time_precision='s')

if __name__ == "__main__":
    main()
