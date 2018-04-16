import logging
import timebuddy
import os
from pymongo import UpdateOne
from copy import deepcopy

monconf_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'
logging.basicConfig(filename='/hpc/local/CentOS7/dhl_ec/software/ctapp/logs/gtex_annot.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger(__name__)

taskid = os.environ['SGE_TASK_ID']
taskid = 1

tissuemap = {1: 'Artery_Aorta',
               2: 'Artery_Coronary',
               3: 'Artery_Tibial',
               4: 'Cells_EBV-transformed_lymphocytes',
               5: 'Cells_Transformed_fibroblasts',
               6: 'Heart_Atrial_Appendage',
               7: 'Heart_Left_Ventricle',
               8: 'Heart_Left_Ventricle',
               9: 'Whole_Blood'}


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


def init_GTEx_annot():
    """
    Update GWAS documents with tissue specific expression data
    First: pull tissue data van GTEx db per tissue (chunked)
    Second: build queries with filters on SNP id (if available)
    or variant id.
    Third: update SNP records with tissue data
    Fourth: Update Gene records with tissue data
    :param records:
    :param col:
    :return:
    """
    cli = get_mongo_client()
    gtex_col = init_db(cli, 'GTEx', 'v7')
    records = get_records(gtex_col)
    colname = ''.join(timebuddy.date().split('-'))
    # get SNP collection object
    gwas_col = init_db(cli, 'MetaGWAS', 'snps{}'.format(colname))
    # get gene collection object
    gene_col = init_db(cli, 'MetaGWAS', 'genes{}'.format(colname))
    update_snp_records(records, gwas_col)
    update_gene_records(records, gene_col)


# write pipeline to database
def bulk_insert(collection, pipeline, name):
    try:
        logger.info('len pipe {}= {}'.format(name, len(pipeline)))
        if len(pipeline) > 1:
            collection.insert_many(pipeline)
    except Exception as e:
        logger.error(e)


def get_records(collection):
    tissue = tissuemap.get(taskid)
    try:
        records = list(collection.find({'tissue': tissue}))
        return records
    except Exception as e:
        logger.error(e)


def update_snp_records(records, col):
    pipeline = list()
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    try:
        for doc in records:
            # make a copy of doc to serve as a update query
            # remove multiple ids from update query
            update = deepcopy(doc)
            update.pop('variant_id')
            update.pop('_id')
            if 'uniprot_id' in doc:
                update.pop('uniprot_id')
            if 'rs_id' in doc:
                update.pop('rs_id')
                snp_update = {'$addToSet': {'tissues': update}}
                snp = doc.get('rs_id')
                query = UpdateOne({'SNP': snp}, snp_update)
            else:
                # else: build queries to match a variant id in the format:chr1:8329030:A:G
                snp_list = list()
                var = doc.get('variant_id')
                vrsplt = var.split('_')
                chr, bp, a1, a2 = vrsplt[0], vrsplt[1], vrsplt[2], vrsplt[3]
                com_a1 = ''.join([complements.get(nuc) for nuc in a1])
                com_a2 = ''.join([complements.get(nuc) for nuc in a2])
                snp_list.append('chr{}:{}:{}:{}'.format(chr, bp, a1, a2))
                # snp_list.append('chr{}:{}:{}:{}'.format(chr, bp, a2, a1))
                snp_list.append('chr{}:{}:{}:{}'.format(chr, bp, com_a1, com_a2))
                # snp_list.append('chr{}:{}:{}:{}'.format(chr, bp, com_a2, com_a1))

                snp_update = {'$addToSet': {'tissues': update}}
                query = UpdateOne({'$in': snp_list}, update)
            pipeline.append(query)
        bulk_insert(col, pipeline, 'SNP GTEx update')
    except Exception as e:
        logger.error(e)


def update_gene_records(records, col):
    pipeline = list()
    try:
        for doc in records:
            update = deepcopy(doc)
            update.pop('variant_id')
            update.pop('_id')
            if 'snp_id' in doc:
                update.pop('snp_id')
            if 'uniprot_id' in doc:
                update.pop('uniprot_id')
                query = {'uniprot_id': doc['uniprot_id']}

                gene_update = {'$addToSet': {'tissues': update}}
                pipeline.append(UpdateOne(query, update))
        bulk_insert(col, pipeline, 'gene GTEx update')
    except Exception as e:
        logger.error(e)
