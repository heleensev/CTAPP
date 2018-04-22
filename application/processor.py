#!/bin/bash
#$ -S /hpc/local/CentOS7/dhl_ec/software/sevpy/1/bin/python
#$ -cwd
#$ -l h_rt=00:05:00
#$ -l h_vmem=5G
#$ -o ~/logs/deploy/
#$ -e ~/logs/deploy/
#$ -N processor

from glob import glob
import os
import subprocess
import itertools
import sys

import fire

base = '/hpc/local/CentOS7/dhl_ec/software/ctapp'
chunker = base+'/application/map_chunker.py'
meta_analysis = 'GWASParser/meta_analysis.py'
gene_mapper = 'gene_mapper.py'
gene_mongo = 'MongoParsenStore/gene_db.py'
geneid_mongo = 'MongoParsenStore/geneID_db.py'
drugbank_mongo = 'MongoParsenStore/drugbank_db.py'
gtex_mongo = 'MongoParsenStore/gtex_db.py'
call_plotter = base+'/application/processor.py'
qq_plotter = ''
mht_plotter = ''

sys.path.append(base)

from application.JobHandler import job_computer
from application.JobHandler import mongo_handler
from application.MetaReader import reader
from application.MetaReader import writer
from application.GWASParser import column_classifier
from application.GWASParser import val_checker
from application.GWASParser import liftover
from application.GWASParser import ref_checker
from application.GWASParser import GWASio
from application.GWASParser import error_logger
from application.MetaReader import config
from application.Annotator import GWAS_db

logpath = os.path.join('/hpc/local/CentOS7/dhl_ec/software/ctapp/', 'logs/processor.log')
monconf_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'
plots_datadir = '/hpc/local/CentOS7/dhl_ec/software/ctapp/data/plots/input'
plots_out = '/hpc/local/CentOS7/dhl_ec/software/ctapp/data/plots/output'

class Processor:

    def __init__(self, conf):
        conf = os.path.abspath(conf)
        self.conf = conf


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

    @staticmethod
    def __get_collection(client, dbname, colname):
        db = client[dbname]
        collection = db[colname]
        return collection

    def preparse(self):

        """
        Method for parsing the GWAS summary csv's and classifying the columns
        before parsing, the config file is validated
        This method is called directly from bash
        :return:
        """
        # validate both the study and the parameters config
        config.ConfigValidator(self.conf, True)

        # build 'study' object with metadata as attributes
        # return type is a list of these study objects
        studies = reader.meta_data(self.conf, 'study')
        params = reader.meta_params(self.conf)

        results = list()
        for study in studies:
            results.append(column_classifier.init_classifier(study))
        # compose new config file with the study objects
        writer.build_config(results, params)

        # log any warnings or errors for user to examine
        error_logger.init_usr_check(results)

    def pipe_execute(self):
        """
        method to launch the jobs on multiple compute nodes via map_chunker.py
        to chunk the individual GWAS study datasets before processing
        :return:
        """
        # validate both the study and the parameters config
        config.ConfigValidator(self.conf, False)
        studies = reader.meta_data(path=self.conf, ob='study')
        params = reader.meta_params(path=self.conf)

        # compute number of jobs to send out, based on config specifications and file size
        mongo_handle = mongo_handler.MongoHandler()
        mongo_handle.launch_mongo(params.mongo, params.runtime, params.mon_vmem)

        # drop collections before rerun
        # cli = self.__get_mongo_client()
        # colname = ''.join(timebuddy.date().split('-'))
        # meta_col = self.__get_collection(cli, 'MetaGWAS', 'meta{}'.format(colname))
        # genes_col = self.__get_collection(cli, 'GWASgenes', 'genes{}'.format(colname))
        # meta_col.drop()
        # genes_col.drop()
        # compute all the possible combinations
        combinations = list()
        for i in range(2, len(studies) + 1):
            comb = list(itertools.combinations(studies, i))
            combinations.extend([list(c) for c in comb])

        for n, study in enumerate(studies):
            job = job_computer.ComputeJobParam(study.path, study.header_indices, study.skip,
                                               params.chunksize, params.subchunk, params.method)
            num_jobs = job.compute()

            # send array job to queue for every computed chunksize in num_job to initiate the pipeline
            subprocess.call('qsub -t 1-{0} -N pre_{1} -l h_rt={2} -l h_vmem={3} {4} {5} '
                            '--config={6} --index={7}'
                            .format(num_jobs, study.studyID, params.runtime, params.vmem, chunker,
                                    params.method, self.conf, [n]), shell=True)
            print('initiating preprocessing for {}'.format(study.studyID))
            # call gene_analysis jobs, waits on preprocessing before proceeding
            subprocess.call('qsub -hold_jid pre_{0} -N MAGMA_{0} -l h_rt=00:10:00 -l h_vmem=15G '
                            '{1} gene_analysis {0} {2} {3}'
                            .format(study.studyID, gene_mapper, study.n_studies, study.ethnicity), shell=True)
        if len(studies) > 1:
            # perform meta analysis on GWAS information already written to the database and update
            subprocess.call('qsub -hold_jid pre_* -l h_rt=00:20:00 -l h_vmem=15G -t 1-24 {0} init_meta_analysis'
                            .format(meta_analysis), shell=True)
            subprocess.call('qsub -hold_jid meta_analysis -l h_rt=00:20:00 -l h_vmem=15G 1-24 -t 1-{} {} plots {} {}'
                            .format((len(studies)+len(combinations)), call_plotter, studies, combinations))
            # # call MAGMA meta_analysis, (per chromosome chunk) wait on preprocessing operation
            subprocess.call('{0} meta_analysis_chunker'.format(gene_mapper), shell=True)

        # do GTEx and Drugdb annotations on GWASdb

        # mongo will be shut down after all the jobs complete
        mongo_handle.stop_mongo(params.mongo, ['meta_analysis, MAGMA_meta*'])

    @staticmethod
    def GWAS_preprocess(data, conf, index, client):

        """
        Method for parsing and pre-processing the GWAS data sets
        (to be executed after the preparse method)
        This method is called via map_chunker.py
        :param data: pandas DataFrame object (chunk of total data set)
        :param index: index for specific document in config file
        :param client: Mongo client for database connection
        :return:
        """
        # validate both the study and the parameters config
        config.ConfigValidator(conf, False)
        # build 'Study' object with metadata as attributes
        study = reader.meta_data(path=conf, ob='study', ix=index[0])
        # check all values in the data set
        data = val_checker.check_correct(data, study)
        # perform liftover on rs ids
        data = liftover.init_liftover(data, client)
        # reference check dataset on 1000 genome data
        data = ref_checker.init_ref_check(data, study, client)
        # write GWAS summary data to the database
        GWAS_db.init_GWAS_db(data, study, client)
        # write to a csv file, for later processing with MAGMA
        GWASio.GWAS_out(data, study.studyID)
        client.close()

    @staticmethod
    def plots(no_studies, combinations):
        # get task ID for index (parallel processing)
        ix = os.environ["SGE_TASK_ID"]
        # get all filenames of plot format files
        qq_files = sorted(glob(os.path.join(plots_datadir, 'qq_*.csv')))
        mht_files = sorted(glob(os.path.join(plots_datadir, 'mht_*.csv')))

        qqfile = qq_files[ix-1]
        name = qqfile.split('qq_')[1].split('.csv')[0]
        # if "_" in the name, Z score of multiple studies was computed, so stat = Z
        stat = 'Z' if '_' in name else 'P'
        # call qqplot.R as subprocess
        subprocess.call('{} -p {} -r {} -o {} -s {} -f PNG '
                        .format(qq_plotter, plots_datadir, qqfile, plots_out, stat), shell=True)
        # if index is smaller or equal to amount of mht files, make both plots else just the QQ plot
        # (for the combination of studies, there is no QQ plot to make)
        if len(mht_files) >= ix:
            mht_file = mht_files[ix-1]
            name = mht_file.split('mht')[1].split('.csv')[0]
            # call manhattan.R as subprocess
            subprocess.call('{} -p {} -r {} -o {} -c FULL -f PNG -t QQ {} '
                            .format(mht_plotter, plots_datadir, mht_file, plots_out, name), shell=True)


    def validate(self):
        """
        method for just the validation of edited config file
        :return:
        """
        config.ConfigValidator(self.conf, True)
        
    def resource_create(self, indices):

        """
        method for creating db collections from the resources: liftover and 1000genome
        :return:
        """
        data = reader.meta_data(self.conf, 'resource', 0)
        params = reader.meta_params(self.conf)

        mongo_handle = mongo_handler.MongoHandler()
        mongo_handle.launch_mongo(params.mongo, params.runtime, params.mon_vmem)

        cli = self.__get_mongo_client()
        ref_col = self.__get_collection(cli, '1000g', 'phase3_2')

        ref_col.drop()

        job = job_computer.ComputeJobParam(data.path, data.header_indices, data.skip,
                                           params.chunksize, params.subchunk, params.method)
        num_jobs = job.compute()

        subprocess.call('qsub -t 1-{} -N create -l h_rt={} -l h_vmem={} {} {} --config={} --index="{}"'
                        .format(num_jobs, params.runtime, params.vmem,
                                chunker, params.method, self.conf, indices), shell=True)

        # mongo will be shut down after all the jobs complete
        mongo_handle.stop_mongo(params.mongo, ['create'])

    def resource_create_bulk(self, database, jobs=None):
        """
        method for the creation of the resources: drugbank and geneID database
        :return:
        """
        params = reader.meta_params(self.conf)

        mongo_handle = mongo_handler.MongoHandler()
        mongo_handle.launch_mongo(True, '00:45:00', '30G')

        creators = {'geneID': geneid_mongo, 'drugbank': drugbank_mongo,
                    'gtex': gtex_mongo}
        db_creator = creators.get(database)
        subprocess.call('qsub -t 1-{} -N gtex -l h_rt={} -l h_vmem={} {}'
                        .format(str(jobs), params.runtime, params.vmem, db_creator))

        mongo_handle.stop_mongo(True, [''])


if __name__ == "__main__":
    fire.Fire(Processor)