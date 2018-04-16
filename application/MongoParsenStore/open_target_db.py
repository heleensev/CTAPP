#!/bin/bash
#$ -S /hpc/local/CentOS7/dhl_ec/software/sevpy/1/bin/python
#$ -cwd
#$ -l h_rt=00:40:00
#$ -l h_vmem=20G
#$ -o ~/logs/deploy/
#$ -e ~/logs/deploy/
#$ -N opentarget

import simplejson as json
from pymongo import ASCENDING, errors
import logging

logging.basicConfig(filename='/hpc/local/CentOS7/dhl_ec/software/ctapp/logs/opentarget.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger(__name__)

monconf_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'
opentarget_file = '/hpc/local/CentOS7/dhl_ec/software/ctapp/data/association_data_heart.json'

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


def bulk_exec(col, pipeline):
    try:
        col.insert_many(pipeline)
    except errors.BulkWriteError as bwe:
        print(bwe.details)


def init_target_db():

    cli = get_mongo_client()
    collection = init_db(cli, 'targetdb', 'opentarget')
    collection.create_index([('geneid', ASCENDING)])
    collection.create_index([('target.uniprot_id', ASCENDING)])
    js_data = open(opentarget_file, 'r')

    pipeline = parse_json(js_data)
    bulk_exec(collection, pipeline)


def collect_ids():
    # get mongo client
    client = get_mongo_client()
    # get geneID collection, for annotation of uniprot ids
    geneID_col = init_db(client, 'IDmapping', 'hgnc_ids')
    try:
        id_annot = list(geneID_col.find({'ensembl': {'$exists': True}}))
        annot_dict = dict()
        for gene_doc in id_annot:
            if gene_doc['ensembl']:
                annot_dict[gene_doc['ensembl']] = gene_doc['uniprot']
        return annot_dict

    except errors.BulkWriteError as be:
        logger.info(be)
    except Exception as e:
        logger.error(e)


def parse_json(js_data):

    uniprot_ids = collect_ids()
    pipeline = list()

    print('parsing json')
    line = True
    while line:
        if line:
            try:
                line = js_data.readline()
                doc = json.loads(line)
                db_doc = dict()
                db_doc['target_id'] = doc['id']
                db_doc['is_direct'] = doc['is_direct']
                db_doc['count'] = doc['evidence_count']['total']
                db_doc['evidence'] = doc['association_score']['datatypes']
                db_doc['score'] = doc['association_score']['overall']
                db_doc['disease'] = doc['disease']['efo_info']['label']
                db_doc['therapeutic_areas'] = doc['disease']['efo_info']['therapeutic_area']['labels']
                target_doc = doc['target']
                uniprot_id = uniprot_ids.get(target_doc['id'])
                if uniprot_id:
                    target_doc['uniprot_id'] = uniprot_id
                db_doc['target'] = target_doc

                pipeline.append(db_doc)

            except json.JSONDecodeError as e:
                "continue to next line"
                print(line)
                print(e)
            except IndexError as e:
                print(e)

    return pipeline


init_target_db()


