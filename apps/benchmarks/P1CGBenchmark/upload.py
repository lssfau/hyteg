import os
import time
import json
from influxdb import InfluxDBClient
from git import Repo


def createFields(data,prefix):
    if not data:
        return []

    json_body = []
    for measure in data:
        if not data[measure]['childs']:
            leaf = True
        else:
            leaf = False
        json_body += [
            {
                'measurement': 'P1_CG_Benchmark',
                'tags': {
                    'host': os.uname()[1],
                    'leaf': leaf,
                },
                'time': int(time.time()),
                'fields':  {prefix + measure: data[measure]['total']}

            }
        ]
        json_body = json_body + createFields(data[measure]['childs'], (prefix + measure + '.'))

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

    with open("P1CGBenchmarkOutput.json") as f:
        data = json.load(f)

    # print(data)
    json_body = createFields(data, '')

    # print(json_body)
    client.write_points(json_body, time_precision='s')


if __name__ == "__main__":
    main()
