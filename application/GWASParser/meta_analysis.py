#!/bin/bash
#$ -S /hpc/local/CentOS7/dhl_ec/software/sevpy/1/bin/python
#$ -cwd
#$ -l h_rt=00:05:00
#$ -l h_vmem=10G
#$ -o ~/logs/deploy/
#$ -e ~/logs/deploy/
#$ -N meta_analysis

from pymongo import UpdateOne
from pymongo import DESCENDING
from pymongo.errors import BulkWriteError
from scipy import stats
from numpy import median
import timebuddy
import platform
import logging
import math
import fire
import re
import os

logpath = '/hpc/local/CentOS7/dhl_ec/software/ctapp/logs'
path = os.path.join(logpath, 'meta_analysis.log')

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

handler = logging.FileHandler(path)
handler.setLevel(logging.INFO)

formatter = logging.Formatter("%(threadName)s:%(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.info('starting meta analysis at {}'.format(timebuddy.time()))

monconf_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'


class MetaAnalysis:

    @staticmethod
    def __get_taskid():
        try:
            taskid = os.environ["SGE_TASK_ID"]
            return taskid
        except AttributeError:
            return False

    @staticmethod
    def __get_mongo_client():
        from pymongo import MongoClient

        def get_host_info():
            with open(monconf_path, 'r') as mon_conf:
                lines = mon_conf.readlines()
            return lines[1].split('=')[1].strip()

        monip = get_host_info()
        client = MongoClient(host=monip, connect=False)
        return client

    def init_meta_analysis(self):
        """
        initiate meta analysis of previously inserted SNPs in the db
        :return:
        """
        client = self.__get_mongo_client()
        chrm = int(self.__get_taskid())

        # collection name is the date of runtime like: 20180130
        colname = ''.join(timebuddy.date().split('-'))
        col = self.init_db(client, 'MetaGWAS', 'meta{}'.format(colname))
        # process SNP per chromosome info in additional chunks of 20
        for no in range(0, 20):
            dbchunk = self.fetch_documents(col, chrm, no)
            # if db output is not empty, ie there are SNPs for this
            # chromosome number
            if dbchunk:
                bulk_pipe = self.update_docs(dbchunk)
                self.bulk_exec(col, bulk_pipe)

    def init_db(self, client, dbname, colname):
        db = client[dbname]
        collection = db[colname]
        return collection

    def bulk_exec(self, col, bulk_obj):
        try:
            col.bulk_write(bulk_obj)
        except BulkWriteError as bwe: 
            logger.info(bwe.details)

    def fetch_documents(self, col, chrm, chnk_no):
        """
        db.products.find({Price:{$gt:290,$lt:400}})
        :param col:
        :param chrm:
        :param chnk_no:
        :return:
        """
        try:
            # db for bp_max should be ref_db not GWASdb
            bp_cur = list(col.find({"CHR": chrm}).sort([('BP', DESCENDING)]).limit(1))
            if len(bp_cur) > 0:
                bp_max = bp_cur[0]['BP']
            else:
                return False
            # get the chunk corresponding to the chnk_no parameter
            under = int((bp_max/20) * chnk_no)
            upper = int(under+(bp_max/20))+1
            chunk = list(col.find({"CHR": chrm, "BP": {'$gt': under, '$lt': upper}}))
            if len(chunk) > 0:
                return chunk
            return False
        except Exception as e:
            logger.error("Error occured during fetch_documents: {}".format(e))


    def update_docs(self, dbchunk):
        """
                ### inverse variance weighted z-score
            $signed_beta[$study] = $sign * $beta[$study];
            $weight[$study] = 1 / ( $se[$study] * $se[$study] );
            my $weighted_beta = $signed_beta[$study] * $weight[$study];
            $total_weighted_beta += $weighted_beta;
            $total_weight += $weight[$study];
            $total_weight_squared += $weight[$study] * $weight[$study];

            ### sample-size weighted z-score
            my $z_weight = sqrt( $sample_size_eff[$study] / $n_eff );
            my $z = ( $signed_beta[$study] / $se[$study] );
            $z_sqrtn += ($z * $z_weight);

            ### sample-size weighted allele frequency
            my $af_weight = $sample_size_eff[$study] / $n_eff;

            doc sample
           {
              "SNP": 278478348,
              "CHR": 1,
              "BP": 72847823,
              "REF": "AT",
              "VAR":
                      [
                        {
                          "ALT": "TC",
                          "STUDY":
                                  [
                                    {
                                      "ID":  "GWAS01",
                                      "FRQ": 0.12,
                                      "SE": 0.34,
                                      "BETA": 0.34
                                    },
                                    {
                                      "ID":  "GWAS02",
                                      "FRQ": 0.13,
                                      "SE": 0.26,
                                      "BETA": 0.24
                                    }
                                  ]
                         "META":
                                {
                                    "Z": "1.34",
                                    "swZ": "1.36",
                                    "P": "1.34",
                                    "swP": "0.17",
                                    "BETA": "0.28",
                                    "SE": "0.21"
                                }
                      ]
          }

        :param dbchunk:
        :return:
        """

        def compute_z():

            # get all squared standard deviations and compute sum

            SEs = [doc['SE'] for doc in var['STUDY']]
            sqSEs = [se**2 for se in SEs]

            # get all the betas in one list and compute sum
            betas = [doc['BETA'] for doc in var['STUDY']]
            wbetas = [(beta / sqSE) for beta, sqSE in zip(betas, sqSEs)]
            wbetasum = sum(wbetas)

            # get all inverted squared standard deviations
            invsqSE = [(1 / sqSE) for sqSE in sqSEs]

            invsqSEsum = sum(invsqSE)

            # compute weighted beta divided by squared standard error
            # and the pooled inverted standard error
            pooled_beta = wbetasum / invsqSEsum
            pooled_se = math.sqrt(1 / invsqSEsum)

            # compute z with the pooled inverse variance-weighted
            # beta coefficient and standard error
            meta_z = pooled_beta / pooled_se

            n_all = [doc['N'] for doc in var['STUDY']]

            n_sum = sum(n_all)
            n_weights = [math.sqrt(n/n_sum) for n in n_all]
            beta_div_SEs = [(beta/SE) for beta, SE in zip(betas, SEs)]
            sample_z = sum([bs * nw for bs, nw in zip(beta_div_SEs, n_weights)])

            # p values for both the standard Z score and the sample size weighted Z score
            p_val = stats.norm.sf(abs(meta_z)) * 2
            sample_p_val = stats.norm.sf(abs(sample_z)) * 2

            # genomic inflation factor, calculated from the sample weighted z score
            chisq = sample_z**2
            #median(chisq) / stats.chisquare(0.5, ddof=1)

            return meta_z, sample_z, pooled_beta, pooled_se, p_val, sample_p_val

        bulk_update = list()
        for doc in dbchunk:
            qfilter = dict()
            qfilter['SNP'] = doc['SNP']
            qfilter['BP'] = doc['BP']

            vars = list(doc['VAR'])
            for var in vars:
                if 'STUDY' in var:
                    if len(var['STUDY']) > 1:
                        z, sample_z, w_beta, w_se, pval, s_pval = compute_z()
                        # add the returned meta analysis values to the doc
                        # format all the numeric values on two decimals
                        var['META'] = {}
                        var['META']['Z'] = '{0:.2f}'.format(z)
                        var['META']['swZ'] = '{0:.2f}'.format(sample_z)
                        var['META']['P'] = '{0:.2f}'.format(z)
                        var['META']['swP'] = '{0:.2f}'.format(s_pval)
                        var['META']['BETA'] = '{0:.2f}'.format(w_beta)
                        var['META']['SE'] = '{0:.2f}'.format(w_se)
                    else:
                        var['META'] = None
                else:
                    vars.remove(var)

            var_doc = {"$set": {'VAR': vars}}
            # add to bulk pipeline
            bulk_update.append(UpdateOne(qfilter, var_doc, upsert=True))

        return bulk_update


if __name__ == "__main__":
    fire.Fire(MetaAnalysis)

    # chunk= [{"BP" : 1036959, "CHR" : 1, "REF" : "C", "SNP" : 11579015, "VAR" : [ { "ALT" : "T", "STUDY" : [ { "ID" : "GWAS01", "FRQ" : 0.10640000000000005, "BETA" : 2.282782465697866, "SE" : 0.165, "N" : 200 }, { "ID" : "GWAS02", "FRQ" : 0.10640000000000005, "BETA" : 2.282782465697866, "SE" : 0.165, "N" : 160 }, { "ID" : "GWAS02", "FRQ" : 0.10640000000000005, "BETA" : 2.282782465697866, "SE" : 0.165, "N" : 160 }, { "ID" : "GWAS01", "FRQ" : 0.10640000000000005, "BETA" : 2.282782465697866, "SE" : 0.165, "N" : 200 } ] }, { "ALT" : "T" } ] }]
    # #
    # m = MetaAnalysis()
    # m.update_docs(chunk)
    """
    { "_id" : ObjectId("5a8dfce7e9d32bb85055d8c5"), "BP" : 1060608, "CHR" : 1, "REF" : "A", "SNP" : 17160824, "VAR" : [ { "ALT" : "G", "STUDY" : [ { "ID" : "GWAS02", "FRQ" : 0.1154, "BETA" : 2.875286120478124, "SE" : 0.1694, "N" : 160 }, { "ID" : "GWAS01", "FRQ" : 0.1154, "BETA" : 2.875286120478124, "SE" : 0.1694, "N" : 200 }, { "ID" : "GWAS02", "FRQ" : 0.1154, "BETA" : 2.875286120478124, "SE" : 0.1694, "N" : 160 }, { "ID" : "GWAS01", "FRQ" : 0.1154, "BETA" : 2.875286120478124, "SE" : 0.1694, "N" : 200 } ], "META" : { "Z" : "33.95", "swZ" : "33.89", "P" : "33.95", "swP" : "0.00", "BETA" : "2.88", "SE" : "0.08" } }, { "ALT" : "G", "META" : null } ] }
    """