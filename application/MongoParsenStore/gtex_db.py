#!/bin/bash
#$ -S /hpc/local/CentOS7/dhl_ec/software/sevpy/1/bin/python
#$ -cwd
#$ -l h_rt=00:15:00
#$ -l h_vmem=15G
#$ -o ~/logs/deploy/
#$ -e ~/logs/deploy/
#$ -N gtexdb

from pymongo import ASCENDING, errors
from sys import argv
import pandas as pd
import logging
import os

logging.basicConfig(filename='/hpc/local/CentOS7/dhl_ec/software/ctapp/logs/gtex.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger(__name__)

taskid = os.environ['SGE_TASK_ID']

gtexdir = '/hpc/local/CentOS7/dhl_ec/software/GTEx/GTEx_Analysis_v7_eQTL/'
idmap = '/hpc/local/CentOS7/dhl_ec/software/GTEx/gtex_variantid_mapping.txt'
monconf_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'

filemap = {1: 'Artery_Aorta.v7.signif_variant_gene_pairs.txt.gz',
           2: 'Artery_Coronary.v7.signif_variant_gene_pairs.txt.gz',
           3: 'Artery_Tibial.v7.signif_variant_gene_pairs.txt.gz',
           4: 'Cells_EBV-transformed_lymphocytes.v7.signif_variant_gene_pairs.txt.gz',
           5: 'Cells_Transformed_fibroblasts.v7.signif_variant_gene_pairs.txt.gz',
           6: 'Heart_Atrial_Appendage.v7.signif_variant_gene_pairs.txt.gz',
           7: 'Heart_Left_Ventricle.v7.signif_variant_gene_pairs.txt.gz',
           8: 'Heart_Left_Ventricle.v7.signif_variant_gene_pairs.txt.gz',
           9: 'Whole_Blood.v7.signif_variant_gene_pairs.txt.gz'}


def get_mongo_client():
    from pymongo import MongoClient

    def get_host_info():
        with open(monconf_path, 'r') as mon_conf:
            lines = mon_conf.readlines()
        return lines[1].split('=')[1].strip()

    monip = get_host_info()
    client = MongoClient(host=monip, connect=False)
    return client


def init_operations():
    datanum = int(taskid)
    filepath = gtexdir+filemap.get(datanum)
    client = get_mongo_client()
    collection = init_db(client, dbname='GTEx', colname='v7')
    # create index for database on "SNP"
    collection.create_index([("rs_id", ASCENDING)])
    collection.create_index([("variant_id", ASCENDING)])
    collection.create_index([('gene_id', ASCENDING)])
    idmap = read_idmap()

    pipeline = process_rows(filepath, idmap)
    bulk_insert(collection, pipeline, 'v7')


def init_db(client, dbname, colname):
    db = client[dbname]
    collection = db[colname]
    return collection


# write pipeline to database
def bulk_insert(collection, pipeline, name):
    try:
        logger.info('len pipe {}= {}'.format(name, len(pipeline)))
        if len(pipeline) > 1:
            collection.insert_many(pipeline)
    except Exception as e:
        logger.error(e)


def read_idmap():
    table = open(idmap, 'r')
    table_doc = {}
    line = True
    while line:
        line = table.readline()
        # check if line not empty
        if line:
            line = line.split()
            var = line[0]
            snp = line[1]
            # add var id to snp id mapping to dict
            table_doc.update({var: snp})

    return table_doc


def collect_ids(genes_list):
    # get mongo client
    client = get_mongo_client()
    # get geneID collection, for annotation of uniprot ids
    geneID_col = init_db(client, 'IDmapping', 'hgnc_ids')
    try:
        print(genes_list[:10])
        id_annot = list(geneID_col.find({'ensembl': {'$in': genes_list}}))
        print(len(id_annot))
        annot_dict = dict()
        for gene_doc in id_annot:
            if gene_doc['ensembl']:
                annot_dict[gene_doc['ensembl']] = gene_doc['uniprot']
        return annot_dict

    except errors.BulkWriteError as be:
        logger.info(be)
    except Exception as e:
        logger.error(e)


def process_rows(data, idmap):

    filename = os.path.basename(data)
    tissue = filename.split('.')[0]
    gtex_df = pd.read_csv(data, header=0, sep='\t', compression='gzip')
    # truncate gene ids (remove .version_number)
    gtex_df[['gene_id']] = gtex_df['gene_id'].dropna().apply(lambda x: x.split('.')[0])
    # convert gene_id column to a list
    all_genes = list(set(gtex_df['gene_id']))
    # get uniprot ids mapping to the ensembl ids in chunks
    uniprotids = collect_ids(all_genes)
    if not uniprotids:
        logger.error('No uniprot ids found, check query and logs')

    # get subset with relevant headers
    gtex_df = gtex_df[['variant_id', 'gene_id', 'pval_nominal', 'slope', 'slope_se']]
    float_cols = ["pval_nominal", "slope", "slope_se"]
    gtex_df[float_cols] = gtex_df[float_cols].astype(float)

    size = gtex_df.shape[0]
    gtex_df['tissue'] = pd.Series([tissue for _ in range(size)])

    gtex_docs = gtex_df.to_dict(orient='records')

    # update gtex records with a rs id and uniprot id, if exists
    for doc in gtex_docs:
        rsid = idmap.get(doc['variant_id'])
        if rsid:
            doc.update({'rs_id': rsid})
        uniprot_id = uniprotids.get(doc['gene_id'])
        if uniprot_id:
            doc.update({'uniprot_id': uniprot_id})

    return gtex_docs


if __name__ == '__main__':
    init_operations()
