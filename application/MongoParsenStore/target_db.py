import simplejson as json
from pymongo import ASCENDING
from glob import glob
from pprint import pprint
import os

monconf_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'

opentarget_dir = '../data/opentarget/'

def get_mongo_client():
    from pymongo import MongoClient

    def get_host_info():
        with open(monconf_path, 'r') as mon_conf:
            lines = mon_conf.readlines()
        return lines[1].split('=')[1].strip()

    monip = get_host_info()
    client = MongoClient(host=monip, connect=False)
    return client


def init_db(client, dbname, colname):
    db = client[dbname]
    collection = db[colname]
    return collection


def bulk_exec(col, bulk_obj):
    try:
        col.bulk_write(bulk_obj, ordered=True)
    except Exception as e:
        pass

def init_target_db():
    cli = get_mongo_client()
    collection = init_db(cli, 'opentarget', 'target1')
    collection.create_index([('geneid', ASCENDING)])
    taskid = os.environ['TASK_ID']

    filename = os.path.join(opentarget_dir, '17{}.json'.format(taskid))
    js_data = open(filename, 'r')

    pipeline = parse_json(js_data)
    bulk_exec(collection, pipeline)

def parse_json(js_data):

    pipeline = list()
    line = True
    while line:
        try:
            line = js_data.readline()
            doc = json.loads(line)
            db_doc = dict()
            db_doc['geneid'] = doc['id'].split('-')[0]
            db_doc['is_direct'] = doc['is_direct']
            db_doc['target'] = doc['target']
            db_doc['evidence'] = doc['association_score']['datatypes']
            db_doc['disease'] = doc['disease']
            pipeline.append(db_doc)

        except json.JSONDecodeError as e:
            "continue to next line"
            print(e)
        except IndexError as e:
            print(e)

    return pipeline

#js_data = json.load(open('/home/sevvy/Scripts/deploy/application/db_samples/opentarget.json', 'rb'))
js_data = open('/home/sevvy/Documents/Stage/Datasets/OpenTargets/17/17_0.json', 'r')

parse_json(js_data)


