#!/bin/bash
#$ -S /hpc/local/CentOS7/dhl_ec/software/sevpy/1/bin/python
#$ -cwd
#$ -l h_rt=00:30:00
#$ -l h_vmem=15G
#$ -o ~/logs/deploy/
#$ -e ~/logs/deploy/
#$ -N bulk_insert


import json
from pymongo import ASCENDING, errors
import logging

logging.basicConfig(filename='/hpc/local/CentOS7/dhl_ec/software/ctapp/logs/bulk_insert.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger(__name__)

monconf_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'
drugbank_json = '/hpc/local/CentOS7/dhl_ec/software/ctapp/data/drugbank/drug_bank_5463.json'

def get_mongo_client():
    from pymongo import MongoClient

    def get_host_info():
        with open(monconf_path, 'r') as mon_conf:
            lines = mon_conf.readlines()
        return lines[1].split('=')[1].strip()
    open(drugbank_json, 'r')
    monip = get_host_info()
    client = MongoClient(host=monip)
    return client


def init_db(client, dbname, colname):
    db = client[dbname]
    collection = db[colname]
    return collection


def bulk_exec(col, pipeline):
    try:
        col.insert(pipeline)
    except errors.BulkWriteError as bwe:
        print(bwe.details)


def insert_json():

    cli = get_mongo_client()
    collection = init_db(cli, 'targetdb', 'drugbank2')
    collection.create_index([('drugbank_id', ASCENDING)])
    collection.create_index([('targets.uniprot_id', ASCENDING)])


    with open(drugbank_json, 'r') as db_json:
        drug_doc = json.load(db_json)
        bulk_exec(collection, drug_doc)

if __name__ == "__main__":
    insert_json()