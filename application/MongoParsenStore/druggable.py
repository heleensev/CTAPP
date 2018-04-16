#!/bin/bash
#$ -S /hpc/local/CentOS7/dhl_ec/software/sevpy/1/bin/python
#$ -cwd
#$ -l h_rt=00:30:00
#$ -l h_vmem=15G
#$ -o ~/logs/deploy/
#$ -e ~/logs/deploy/
#$ -N druggable

import logging
import sys
import pandas as pd
import simplejson as json
from pymongo import ASCENDING, errors

base = '/hpc/local/CentOS7/dhl_ec/software/ctapp'
sys.path.append(base)
from application.JobHandler import mongo_handler

logging.basicConfig(filename='/hpc/local/CentOS7/dhl_ec/software/ctapp/logs/druggable.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger(__name__)

monconf_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'
finan_file = '/hpc/local/CentOS7/dhl_ec/software/ctapp/data/aag1166_Table_S1.xlsx'


def launch_mongo():
    mongo_handle = mongo_handler.MongoHandler()
    mongo_handle.launch_mongo(True, '00:01:00', '30G', True)


def init_db(client, dbname, colname):
    db = client[dbname]
    collection = db[colname]
    return collection


def get_mongo_client():
    from pymongo import MongoClient

    def get_host_info():
        with open(monconf_path, 'r') as mon_conf:
            lines = mon_conf.readlines()
        return lines[1].split('=')[1].strip()

    monip = get_host_info()
    client = MongoClient(host=monip, connect=False)
    return client


def init_druggable():
    # # launch database on a node
    # launch_mongo()
    cli = get_mongo_client()
    collection = init_db(cli, 'targetdb', 'finan')
    collection.create_index([('uniprot_id', ASCENDING)])
    fin_csv = pd.read_excel(finan_file, header=0)
    parse_druggable(fin_csv, collection)


def collect_ids():
    # get mongo client
    client = get_mongo_client()
    # get geneID collection, for annotation of uniprot ids
    geneID_col = init_db(client, 'IDmapping', 'hgnc_ids')
    try:
        id_annot = list(geneID_col.find({'uniprot': {'$exists': True}}))
        annot_dict = dict()
        for gene_doc in id_annot:
            if gene_doc['ensembl']:
                annot_dict[gene_doc['ensembl']] = gene_doc['uniprot']
        return annot_dict

    except errors.BulkWriteError as be:
        logger.info(be)
    except Exception as e:
        logger.error(e)


# write pipeline to database
def bulk_insert(collection, pipeline, name):
    try:
        logger.info('len pipe {}= {}'.format(name, len(pipeline)))
        if len(pipeline) > 1:
            collection.insert(pipeline)
    except Exception as e:
        logger.error(e)
        with open('drugbank.json', 'w') as out:
            json.dump(pipeline, out)


def parse_druggable(druggable, collection):

    y_n = {'Y': True,
           'N': False}

    uniprot_ids = collect_ids()

    pipeline = list()
    for i, gene_name, tier, symbol, chrm, start, end, strand, descr, gwas, small, bio, adme in druggable.itertuples():
        druggable_doc = dict()
        druggable_doc['gene_name'] = gene_name
        # check for uniprot id in the list and update if present
        uni = uniprot_ids.get(gene_name)
        if uni:
            druggable_doc['uniprot_id'] = uni
        druggable_doc['description'] = descr
        druggable_doc['gwas_regions'] = int(gwas)
        druggable_doc['small_mol_drug'] = y_n[small]
        druggable_doc['bio_therapeutic'] = y_n[bio]
        druggable_doc['ADME'] = y_n[adme]

        # append to the pipeline
        pipeline.append(druggable_doc)

    bulk_insert(collection, pipeline, 'finan')

if __name__ == "__main__":
    init_druggable()


