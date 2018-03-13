# module parsing
import pandas as pd
import logging, sys
import pymongo
import platform
from pymongo import errors
import warnings
import timebuddy
import os
import re
#warnings.simplefilter(action='ignore', category='FutureWaring')
pd.options.mode.chained_assignment = None

logpath = '/hpc/local/CentOS7/dhl_ec/software/ctapp/logs'

path = os.path.join(logpath, 'liftover_creator.log')

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

handler = logging.FileHandler(path)
handler.setLevel(logging.INFO)

formatter = logging.Formatter("%(threadName)s:%(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.info('starting liftover_creator at: {}'.format(timebuddy.time()))


def init_db(client, dbname, colname):
    db = client[dbname]
    collection = db[colname]
    return collection


def bulk_insert(collection, pipeline, name):
    try:
        # do bulk insert of obsolete rs ids
        logger.info('len pipe {}= {}'.format(name, len(pipeline)))
        if len(pipeline) > 1:
            collection.insert_many(pipeline, ordered=False)
    except pymongo.errors.BulkWriteError as e:
        logger.error(e)
        #nerrs = len(e._OperationFailure__details.get('writeErrors'))
        #logger.error('number of failed writes: {}'.format(nerrs))
    except Exception:
        logger.error(sys.exc_info())
        # inserted = r._InsertManyResult__inserted_ids
        # not_inserted = [doc['_id'] for doc in pipeline if doc not in inserted]
        # logger.error('not inserted: {}'.format(not_inserted))


def bulk_update(collection, pipeline):
    try:
        # do bulk update of obsolete rs ids
        logger.info('len pipe: {}'.format(len(pipeline)))
        collection.update_many(pipeline, ordered=False)
    except errors.BulkWriteError as e:
        logger.error(e)
        nerrs = len(e._OperationFailure__details.get('writeErrors'))
        logger.error('number of failed writes: {}'.format(nerrs))
    except Exception:
        logger.error(sys.exc_info(), sys.exc_info()[-1].tb_lineno)


def bulk_remove(collection, pipeline):
    try:
        collection.remove(pipeline)
    except pymongo.errors.BulkWriteError as e:
        logger.error(e)
        nerrs = len(e._OperationFailure__details.get('writeErrors'))
        logger.error('number of failed writes: {}'.format(nerrs))
    except Exception:
        logger.error(sys.exc_info(), sys.exc_info()[-1].tb_lineno)


def snp_synonyms(df, client):

    """
    function to create a mongo collection containing old rs IDs
    :param df: complete dataframe of RsMergeArch file
    :param client: mongo client
    """
    collection = init_db(client, colname='synonyms', dbname='rsliftover')
    # filter out the duplicate current rs ids and sort them
    df['idx_org'] = df.groupby(['rscur']).rscur.transform('idxmin')
    df = df[df.duplicated(subset=['rscur'], keep='first')]
    df = df.sort_values(by=['rscur'])
    # put all the rs ids that map to the same current rs id in a list
    rsdict = dict()
    last = int(df.iloc[0].rscur)
    rsdict[last] = {int(df.iloc[0].idx_org)}
    for high, cur, idx in df.itertuples():
        high, cur, idx = int(high), int(cur), int(idx)
        if cur != last:
            rsdict[last] = list(rsdict[last])
            last = cur
            rsdict[cur] = {idx}
        rsdict[cur].add(high)
    rsdict[last] = list(rsdict[last])
    # transform the dictionary to a list to perform bulk insert
    doc = [{'rs_cur': k, 'old': v} for k, v in rsdict.items()]

    bulk_insert(collection, doc, 'synonyms')


def snp_liftover(mergearch, rshistory, client):

    """
    function to create a mongo collection with old to new rsID mapping
    :param mergearch: dataframe of complete RsMergeArch data
    :param rshistory: dataframe of complete SNPHistory data
    :param client: mongo client object
    :return:
    """

    collection = init_db(client, colname='mergearch', dbname='rsliftover')

    # gives you all the rows without the rsm ones, but Nans are thrown out
    nonrsm = rshistory[rshistory.comment.str.contains("[^(rsm)]", na=False)]

    # filter rshistory on Nans and create new dataframe
    nans = rshistory[rshistory.comment.isnull()]
    # concat the Nans dataframe with non-rsm dataframe
    rshistory = nonrsm.loc[:, 'rscur'].append(nans.loc[:, 'rscur']).to_frame()

    # get overlap of SNPHistory dataset and RsMergeArch to find retracted SNPs, add bool column
    history = pd.merge(rshistory, mergearch, how='inner', on=['rscur'])
    history.apply(pd.to_numeric, errors='ignore')
    history.loc[:, 'hist'] = pd.Series(True, index=history.index)

    # get the remaining SNPs (non-retracted or nonhistory)

    nonhistory = mergearch[(~mergearch.rscur.isin(history.rscur))]
    nonhistory.apply(pd.to_numeric, errors='ignore')
    nonhistory.loc[:, 'hist'] = pd.Series(False, index=nonhistory.index)

    # convert to dictionary
    hist_pipe = history.to_dict(orient='records')
    nonhist_pipe = nonhistory.to_dict(orient='records')

    bulk_insert(collection, hist_pipe, 'history')
    bulk_insert(collection, nonhist_pipe, 'non_history')
