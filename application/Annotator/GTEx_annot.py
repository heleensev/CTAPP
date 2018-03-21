import logging
import timebuddy
import os
from pymongo import UpdateOne

monconf_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'
logging.basicConfig(filename='../../logs/gtex_annot.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger(__name__)

taskid = os.environ['TaskID']
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

def init_GTEx():

    cli = get_mongo_client()
    gtex_col = init_db(cli, 'GTEx', 'v7')
    records = get_records(gtex_col)
    colname = ''.join(timebuddy.date().split('-'))
    gwas_col = init_db(cli, 'MetaGWAS', 'meta{}'.format(colname))
    update_gwas(records, gwas_col)
    pass

# write pipeline to database
def bulk_update(collection, pipeline, name):
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


def update_gwas(records, col):
    pipeline = list()
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    for doc in records:
        if 'rsid' in records:
            snp = doc.get('rsid')
            query = UpdateOne({'SNP': snp}, {})

        else:
            snp_list = list()
            var = doc.get('variant_id')
            vrsplt = var.split('_')
            "use $in operator!!"
            if len(vrsplt) == 5:
                chr, bp, a1, a2 = vrsplt[0], vrsplt[1], vrsplt[2], vrsplt[3]
                com_a1 = ''.join([complements.get(nuc) for nuc in a1])
                com_a2 = ''.join([complements.get(nuc) for nuc in a2])
                snp_list.append('chr{}:{}:{}:{}'.format(chr, bp, a1, a2))
                snp_list.append('chr{}:{}:{}:{}'.format(chr, bp, a2, a1))
                snp_list.append('chr{}:{}:{}:{}'.format(chr, bp, com_a1, com_a2))
                snp_list.append('chr{}:{}:{}:{}'.format(chr, bp, com_a2, com_a1))


    pass