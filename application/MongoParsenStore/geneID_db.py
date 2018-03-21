#!/bin/bash
#$ -S /hpc/local/CentOS7/dhl_ec/software/sevpy/1/bin/python
#$ -cwd
#$ -l h_rt=00:10:00
#$ -l h_vmem=10G
#$ -o ~/logs/deploy/
#$ -e ~/logs/deploy/
#$ -N geneiddb

import requests
import logging
import pandas as pd
from pymongo import ASCENDING
from io import StringIO
from os import environ

base = '/hpc/local/CentOS7/dhl_ec/software/ctapp/logs/'
logging.basicConfig(filename=base + 'snpid_db.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger(__name__)
monconf_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'


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


def init_geneid_db(chr):
    client = get_mongo_client()
    col = init_db(client, 'IDmapping', 'hgnc_ids2')

    # create indexes for fast querying
    col.create_index([('entrez', ASCENDING)])
    col.create_index([('uniprot', ASCENDING)])
    col.create_index([('ensembl', ASCENDING)])
    col.create_index([('symbol', ASCENDING)])

    # chromosome 23 equals chromosome X
    if chr == '23':
        chr = 'X'

    pipe = parse_records(chr)
    bulk_insert(col, pipe, chr)


# # write pipeline to database
def bulk_insert(collection, pipeline, name):
    try:
        logger.info('len pipe chromosome {}= {}'.format(name, len(pipeline)))
        if len(pipeline) > 1:
            collection.insert_many(pipeline)
    except Exception as e:
        logger.error(e)


def parse_records(chr):
    base_url = "https://www.genenames.org/cgi-bin/download?col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col" \
               "=gd_aliases&col=md_eg_id&col=md_prot_id&col=md_ensembl_id&status=Approved&status=Entry+Withdrawn" \
               "&status_opt=2&chr={}&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit"

    gene_ids = ['symbol', 'name', 'previous', 'synonyms', 'entrez', 'uniprot', 'ensembl']
    r = requests.get(base_url.format(chr))
    content = r.text
    buffer = StringIO(content)
    gene_data = pd.read_csv(buffer, sep='\t', names=gene_ids, skiprows=1, header=None)

    # replace numpy NaNs with empty strings
    gene_data.fillna('', inplace=True)

    # convert dataframe to json type records for database storing
    geneid_doc = gene_data.to_dict(orient='records')

    # convert comma separated fields with lists
    for i, doc in enumerate(geneid_doc):
        syns = doc['synonyms']
        doc['synonyms'] = [s.strip() for s in str(syns).split(',')]
        syns = doc['previous']
        doc['previous'] = [s.strip() for s in str(syns).split(',')]
        geneid_doc[i] = doc

    return geneid_doc

if __name__ == "__main__":
    chr = environ['SGE_TASK_ID']
    init_geneid_db(str(chr))
