import logging
import timebuddy
import os
from pymongo import UpdateOne, errors

monconf_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'
logging.basicConfig(filename='/hpc/local/CentOS7/dhl_ec/software/ctapp/logs/opentarget_annot.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger(__name__)


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
    Update Gene documents with OpenTarget data
    First: Pull
    Second: build queries with filters on SNP id (if available)
    or variant id,
    Third: update SNP records with tissue data
    Fourth: Update Gene records with tissue data
    :param records:
    :param col:
    :return:
    """

    # get mongo client
    client = get_mongo_client()
    colname = ''.join(timebuddy.date().split('-'))
    gene_col = init_db(client, 'MetaGWAS', 'genes{}'.format(colname))

    study_ids = get_study_ids(gene_col)
    for study in study_ids:
        gene_ids = get_study_gene_ids(study, gene_col)
        pipeline = annotate_records(gene_ids)
        bulk_exec(gene_col, pipeline)

def bulk_exec(col, pipeline):
    try:
        col.insert(pipeline)
    except errors.BulkWriteError as bwe:
        print(bwe.details)

def get_study_ids(gene_col):
    try:
        studies = list(gene_col.distinct('STUDY.StudyID'))
        return studies
    except errors.ConnectionFailure as cf:
        logger.error(cf)
    except Exception as e:
        logger.error(e)



def get_study_gene_ids(gene_col, study_id):
    try:
        gene_ids = list(gene_col.find({'STUDY.studyID': study_id, 'uniprot_id': {'$exists': True}}, {'uniprot_id': 1}))
        return gene_ids
    except errors.PyMongoError as pe:
        logger.error(pe)
    except Exception as e:
        logger.error(e)



def match_records(db, uni_map, gene_list, fetch_field, filter_field=False):
    # get mongo client
    client = get_mongo_client()
    # get geneID collection, for annotation of uniprot ids
    geneID_col = init_db(client, 'targetdb', db)
    try:
        query_field = {uni_map: {'$in': gene_list}}
        if filter_field:
            query_field.update(filter_field)
        id_annot = list(geneID_col.find(filter_field, fetch_field))
        annot_dict = dict()
        for gene_doc in id_annot:
            gene_doc.update({'db': db})
            annot_dict[gene_doc['uniprot_id']] = gene_doc

        return annot_dict

    except errors.BulkWriteError as be:
        logger.info(be)
    except Exception as e:
        logger.error(e)


def annotate_records(gene_list):
    pipeline = list()

    # drug target databases
    # corresponding uniprot id mappings
    # corresponding filters or 'where' queries
    # corresponding field to fetch
    annot_dbs = ['drugbank', 'opentarget', 'TTD']
    uniprot_field = ['targets.uniprot_id', 'target.uniprot_id', 'uniprot']
    filter_field = [False, False, {'score': {'$gt': 0.5}}]
    fetch_field = [{"drugbank_id": 1, "name": 1, "status": 1, "pharmacodynamics": 1, "indication": 1},
                   {"target_id": 1, "score": 1, "evidence": 1, "is_direct": 1, "disease":1},
                   {"target_id": 1, "type": 1, "diseases": 1, "function": 1}]

    for db, uni, filtr, fetch in zip(annot_dbs, uniprot_field, filter_field, fetch_field):
        records = match_records(db, uni, gene_list, filtr, fetch)
        if records:

            for uni, annot in records.items():
                filter_query = uni
                update_query = annot
                pipeline.append(UpdateOne(filter_query, update_query))

    return pipeline
