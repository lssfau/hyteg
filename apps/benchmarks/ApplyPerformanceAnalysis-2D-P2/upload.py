import os
import time
import json
from influxdb import InfluxDBClient


def createFields(data):

    if not data:
        return []

    json_body = []

    for measure in data:
        level = int(str(measure[measure.find('level') + 5:measure.find('level') + 7]))
        json_body += [
            {
                'measurement': 'ApplyPerformanceAnalysis-2D-P2-local-bylevel-nochildren2',
                'tags': {
                    'host': os.uname()[1],
                    'level': level,
                },
                'time': int(time.time()),
                'fields':  {measure: data[measure]['total']}
            }
        ]
    return json_body


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

    # repo = Repo(search_parent_directories=True)
    # commit = repo.head.commit

    with open("ApplyPerformanceAnalysis-2D-P2.json") as f:
        data = json.load(f)

    # print(data)
    json_body = createFields(data)

    print(json_body)
    client.write_points(json_body, time_precision='s')


if __name__ == "__main__":
    main()
