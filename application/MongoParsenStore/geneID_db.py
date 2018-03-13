import requests
import logging
import pandas as pd
from pymongo import ASCENDING
from io import StringIO

logging.basicConfig(filename='../logs/snpid_db.log', filemode='a', level=logging.DEBUG)
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


def init_geneid_db():
    client = get_mongo_client()
    col = init_db(client, 'gene_ids', 'hgnc')

    # create indexes for fast querying
    col.create_index([('geneid', ASCENDING)])
    col.create_index([('uniprot')], ASCENDING)
    col.create([('ensembl', ASCENDING)])
    col.create_index(['chr', ASCENDING])

    all_chr = [x for x in range(1, 23)]
    all_chr += ['X', 'Y']
    for chr in all_chr:
        pipe = parse_records(chr)
        bulk_insert(col, pipe, chr)


# # write pipeline to database
def bulk_insert(collection, pipeline, name):
    try:
        logger.info('len pipe {}= {}'.format(name, len(pipeline)))
        if len(pipeline) > 1:
            collection.insert_many(pipeline)
    except Exception as e:
        logger.error(e)


def parse_records(chr):

    base_url = "https://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name" \
               "&col=gd_aliases&col=gd_pub_chrom_map&col=md_eg_id&col=md_prot_id&col=md_ensembl_id" \
               "&status=Approved&status=Entry+Withdrawn&status_opt=2&chr={}&where=" \
               "&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit"

    gene_ids = ['HGNC', 'symbol', 'name', 'synonyms', 'chr', 'entrez', 'geneid', 'uniprot', 'ensembl']
    r = requests.get(base_url.format(chr))
    cont = r.text
    buffer = StringIO(cont)
    gene_data = pd.read_csv(buffer, sep='\t', header=0, names=gene_ids)

    # make sure the entrez ids are stored as integers, not floats
    gene_data[['entrez']] = gene_data['entrez'].fillna('')
    gene_data[gene_ids.remove('entrez')] = gene_data[gene_ids.remove('entrez')].fillna('')
    # convert dataframe to json type records for database storing
    geneid_doc = gene_data.to_dict(orient='records')

    return geneid_doc

if __name__ == "__main__":
    init_geneid_db()