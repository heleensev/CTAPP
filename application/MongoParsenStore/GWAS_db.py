# module parsing
import pandas as pd
import logging, sys
from pymongo import UpdateOne, ASCENDING
from pymongo.errors import BulkWriteError
import platform
from pymongo import errors
import warnings
import timebuddy
import pprint
import os
import re
#warnings.simplefilter(action='ignore', category='FutureWaring')
pd.options.mode.chained_assignment = None

logpath = '/hpc/local/CentOS7/dhl_ec/software/ctapp/logs'

path = os.path.join(logpath, 'gwas_db.log')

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

handler = logging.FileHandler(path)
handler.setLevel(logging.INFO)

formatter = logging.Formatter("%(threadName)s:%(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.info('starting GWAS_db at: {}'.format(timebuddy.time()))


"""
db.col.update({"SNP": 12565280, "CHR": 1, "BP": 711153, "REF": "A"},
{$push: {VAR: {"ALT": "C", "FRQ": 0.0662, "OR": -0.914, "SE": 0.38923}}}, {upsert:true})

db.col.update({"SNP": 12565270, "CHR": 1, "BP": 711153, "REF": "A"},
{$addToSet: {VAR: {"ALT": "C"}}}, {upsert:true})

db.col.update({"SNP": 12565270, "VAR": {$elemMatch: {"ALT":"C"}}},
{$push: {"VAR.$.STUDY": {"FRQ": 0.0662, "OR": -0.914, "SE": 0.38923}}}, {upsert:true})


db.col.update({"SNP": 12565280, "CHR": 1, "BP": 711153, "REF": "A", "VAR": [{'ALT': "C"}]},
{$push: {"VAR$STUDY": {"FRQ": 0.0662, "OR": -0.914, "SE": 0.38923}}}, {upsert:true})
"""

"""
reqs = UpdateOne(filter={"SNP": 12565280, "CHR": 1, "BP": 711153, "REF": "A"},
update={'$push': {'VAR': {"ALT": "CG", "FRQ": 0.4378, "OR": -0.380, "SE": 0.4785}}}, upsert=True)"""


def init_GWAS_db(data, study, client):
    """
    initiator of GWAS_creator: creates collection containing all the meta data of
    the individual GWAS studies

    :param data: dataframe chunk of a GWAS study
    :param study: Study object, with study specific attributes
    :param client: pymongo MongoClient object
    :return:
    """
    rs, no_rs = data
    # collection name is the date of runtime like: 20180130
    colname = ''.join(timebuddy.date().split('-'))
    col = init_db(client, 'MetaGWAS', 'meta{}'.format(colname))
    # create index for collection
    col.create_index([("SNP", ASCENDING)])
    # create json docs and add to a pipeline
    bulk_ob = data_to_json(rs, 'rs', study)
    # push pipeline to database
    bulk_exec(col, bulk_ob)
    # do the same if a dataframe with no SNPs rs IDs exists
    if not no_rs.empty:
        bulk_ob = data_to_json(no_rs, 'no_rs', study)
        bulk_exec(col, bulk_ob)


def init_db(client, dbname, colname):
    db = client[dbname]
    collection = db[colname]
    return collection


def bulk_exec(col, bulk_obj):
    try:
        col.bulk_write(bulk_obj, ordered=True)
    except BulkWriteError as bwe:
        pprint.pprint(bwe.details)
        logger.error(bwe.details)


def data_to_json(data, prefix, study):

    """
    convert the dataframe with GWAS summary data
    to a json format for storing in the database
    :param data: pandas Dataframe
    :param prefix: "rs" or "no_rs"
    :param study: Study object with GWAS specific attributes
    :return:
    """
    ID = study.studyID
    n_studies = study.n_studies
    h = list(data.columns.values)
    # the dataframe may not be in the order of the listed headers below, so indexing is done
    SNP, CHR, BP, EA, OA, FRQ, BETA, SE, P = \
        h.index('SNP'), h.index('CHR'), h.index('BP'), h.index('EA'), h.index('OA'), \
        h.index('FRQ'), h.index('BETA'), h.index('SE'), h.index('P')
    if n_studies:
        n = False
    else:
        n = True
        N = h.index('N')
    bulk = list()
    for row in data.itertuples(index=False):
        # get the right elements from the data frame row on header indices
        snp, chrm, bp, a1, a2, frq, bet, se, p = \
            row[SNP], row[CHR], row[BP], row[EA], row[OA], row[FRQ], row[BETA], row[SE], row[P]

        # first, set the info of this genomic loc if it not exists in the db
        q_filter = dict()
        q_filter['SNP'] = snp
        loc_info = dict()
        loc_info['SNP'] = snp
        loc_info['CHR'] = int(chrm)
        loc_info['BP'] = int(bp)
        loc_info['REF'] = a1
        q_update = {'$setOnInsert': loc_info}
        # append the update to the bulk pipeline
        bulk.append(UpdateOne(q_filter, q_update, upsert=True))

        # second, add the variant (allele) subdocument if it not exists
        q_filter = dict()
        q_filter['SNP'] = snp
        q_filter['VAR.ALT'] = {'$ne': a2}
        q_update = {'$addToSet': {'VAR': {'ALT': a2}}}
        bulk.append((UpdateOne(q_filter, q_update)))

        # third, add the info of particular SNP to the subdocument
        snp_info = dict()
        snp_info['ID'] = ID
        snp_info['FRQ'] = float(frq)
        snp_info['BETA'] = float(bet)
        snp_info['SE'] = float(se)
        # the study size may be a fixed value per study or a unique value per SNP
        if n:
            snp_info['N'] = row[N]
        else:
            snp_info['N'] = int(n_studies)
        # expand query to match the previously upserted alt allele
        # and add the study specific metadata of the SNP
        q_filter = dict()
        q_filter['SNP'] = snp
        q_filter['VAR'] = {'$elemMatch': {'ALT': a2}}
        study_update = {'$push': {"VAR.$.STUDY": snp_info}}

        bulk.append(UpdateOne(q_filter, study_update))
    return bulk

"""
{'$addToSet': {'VAR': {'ALT': 'A'}}}
{'SNP': 'rs2905035', 'CHR': 1, 'BP': 775659, 'REF': 'G'}
{'$push': {'VAR.$.STUDY': {'ID': 'GWAS01', 'FRQ': 0.85970000000000002, 'BETA': 1.6024622986976751, 'SE': 0.2157, 'N': 200}}}
{'SNP': 'rs2905035', 'CHR': 1, 'BP': 775659, 'REF': 'G', 'VAR': {'$elemMatch': {'ALT': 'A'}}}
{'$addToSet': {'VAR': {'ALT': 'A'}}}

db.meta201802115.update({'SNP': 'rs2905035', 'CHR': 1, 'BP': 775659, 'REF': 'G'}, {'$addToSet': {'VAR': {'ALT': 'A'}}}, {upsert:true})

db.meta201802115.update({'SNP': 'rs2905035', 'CHR': 1, 'BP': 775659, 'REF': 'G', 'VAR': {'$elemMatch': {'ALT': 'A'}}}, {'$push': {'VAR.$.STUDY': {'ID': 'GWAS01', 'FRQ': 0.85970000000000002, 'BETA': 1.6024622986976751, 'SE': 0.2157, 'N': 200}}}, {upsert:true})
"""