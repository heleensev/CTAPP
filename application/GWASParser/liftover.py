from pymongo import errors
import simplejson as json
import pandas as pd

from timebuddy import time
import logging ,sys

GWAS = ""
liftover_sum = 0
nohit_sum = 0

logging.basicConfig(filename='../logs/liftover.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger(__name__)


def init_db(client, dbname, colname):
    db = client[dbname]
    collection = db[colname]
    return collection


def init_liftover(data, client):
    rs, no_rs = data
    collection = init_db(client, 'rsliftover', 'merge')
    snp_col = rs['SNP'].tolist()
    id_array = [int(snp.strip('rs')) for snp in snp_col]
    # find matching rs ids in the reference set to be converted
    docs = bulk_find(collection, id_array)
    # convert the rs ids to integers so that they are compatible with the reference data
    rs[['SNP']] = rs['SNP'].apply(lambda x: int(x.strip('rs')))
    # get the SNPs from the database and join results with dataframe
    rs = update_df(docs, rs)
    # convert the ids back to rs ids by prepending "rs"
    rs[['SNP']] = rs['SNP'].apply(lambda x: 'rs{}'.format(str(x)))

    return rs, no_rs


def update_df(docs, data):
    docs = list(docs)
    if len(docs) > 0:
        # dump results from query to a json string (= dictionary-like)
        json_docs = json.dumps(docs)
        lifted = pd.read_json(json_docs, orient='records')

        # merge gwas dataframe with results with sql-like left join
        dfmerge = pd.merge(data[['SNP']], lifted, left_on='SNP', right_on='rshigh', how='left')

        # use boolean masking to replace original rs ids with lifted rs ids fromqst
        # query and keeping the original (unlifted) rs ids
        data.loc[dfmerge['rscur'].notnull(), 'SNP'] = dfmerge['rscur']
    return data


def bulk_find(collection, querylist):
    try:
        results = collection.find({"rshigh": {"$in": querylist}}, {"_id": 0, "rscur": 1, "rshigh": 1})
        return results
    except errors.ConnectionFailure:
        print('ehhh')
