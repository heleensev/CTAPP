#!/bin/bash
#$ -S /hpc/local/CentOS7/dhl_ec/software/sevpy/1/bin/python
#$ -cwd
#$ -l h_rt=00:40:00
#$ -l h_vmem=20G
#$ -o ~/logs/deploy/
#$ -e ~/logs/deploy/
#$ -N ttd_db

import logging
from pymongo import ASCENDING


monconf_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'
logging.basicConfig(filename='/hpc/local/CentOS7/dhl_ec/software/ctapp/logs/ttd_db.log', filemode='a', level=logging.DEBUG)

logger = logging.getLogger(__name__)

full_db = '/hpc/local/CentOS7/dhl_ec/software/ctapp/data/P1-01-TTD_download.txt'


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


def init_ttd_db():
    cli = get_mongo_client()
    col = init_db(cli, 'targetdb', 'TTD')
    col.create_index([('uniprot_id', ASCENDING)])
    col.create_index([('target_id', ASCENDING)])
    parse_full_ttd(col)

# write pipeline to database
def bulk_insert(collection, pipeline, name):
    try:
        logger.info('len pipe {}= {}'.format(name, len(pipeline)))
        if len(pipeline) > 1:
            collection.insert_many(pipeline)
    except Exception as e:
        logger.error(e)


def parse_full_ttd(collection):

    with open(full_db, 'r') as ttd:
        # read header lines
        for i in range(12):
            line = ttd.readline().strip()
        line = True
        target_prev = ''
        target_id, uniprot_id, name, type, function = '', '', '', '', ''
        drugs, diseases, pathways = [], [],  set()
        target_bulk = list()

        def update_doc():
            target_bulk.append(
                {'target_id': target_id,
                 'uniprot': uniprot_id,
                 'target_name': name,
                 'type': type,
                 'function': function,
                 'drugs': drugs,
                 'diseases': diseases,
                 'pathways': list(pathways)})

        while line:
            line = ttd.readline().strip()
            # check if line not empty string
            if line:
                split_line = line.split('\t')
                target_id = split_line[0]
                # check if the target ID is still the same
                # if not: add previous list and dicts to the pipeline
                # and create new list and dict
                if target_id != target_prev:
                    if target_prev:
                        update_doc()
                    target_prev = target_id
                    drugs, diseases, pathways = list(), list(), set()

                if 'UniProt ID' == split_line[1]:
                    uniprot_id = split_line[2]
                elif 'Name' in split_line:
                    name = split_line[2]
                elif 'Type of target' == split_line[1]:
                    type = split_line[2]
                elif 'Function' == split_line[1]:
                    function = split_line[2]
                elif 'Drug(s)' == split_line[1]:
                    drugs.append(split_line[2])
                elif 'Disease' == split_line[1]:
                    diseases.append(split_line[2])
                elif "Pathway" in split_line[1]:
                    pathways.add(split_line[2])
        # append information of last entry to the pipeline
        update_doc()

    bulk_insert(collection, target_bulk, 'TTD')

if __name__== "__main__":
    init_ttd_db()