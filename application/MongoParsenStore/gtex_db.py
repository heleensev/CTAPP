from pymongo import ASCENDING
from sys import argv
import pandas as pd
import logging
import os

logging.basicConfig(filename='../../logs/gtex.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger(__name__)

# taskid = os.environ['TASK_ID']
taskid = 1
gtexdir = '/hpc/local/CentOS7/dhl_ec/software/GTEx/'
idmap = gtexdir+'GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt'
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

    datanum = argv[1]
    filepath = gtexdir+filemap.get(datanum)
    client = get_mongo_client()
    collection = init_db(client, dbname='GTEx', colname='v7')
    # create index for database on "SNP"
    collection.create_index([("rsid", ASCENDING)])
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

    line = table.readline()
    while line:
        line = table.readline().split()
        var = line[2]
        snp = line[6]

        table_doc.update({snp: var})

    return table_doc

def process_rows(data, idmap):

    filename = os.path.basename(data)
    tissue = filename.split('.')[0]
    df = pd.read_csv(data, header=0, sep='\t', compression='gzip')

    # get subset with relevant headers
    df = df[['variant_id', 'gene_id', 'pval_nominal', 'slope', 'slope_se']]
    float_cols = ["pval_nominal", "slope", "slope_se"]
    df[float_cols] = df[float_cols].astype(float)

    size = df.shape[0]
    df['tissue'] = pd.Series([tissue for _ in range(size)])

    gtex_docs = df.to_dict(orient='records')

    # update gtex records with a rsid if present, else a None is added
    for doc in gtex_docs:
        rsid = idmap.get(doc['variant_id'])
        if rsid != 0:
            doc.update({'rsid': rsid})
        else:
            doc.update({'rsid': None})

    return gtex_docs


if __name__ == '__main__':
    init_operations()
