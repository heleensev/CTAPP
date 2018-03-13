#!/bin/bash
#$ -S /hpc/local/CentOS7/dhl_ec/software/sevpy/1/bin/python
#$ -cwd
#$ -o ~/logs/deploy/
#$ -e ~/logs/deploy/
#$ -N mapchunker

import os
import sys
import fire
import yaml
import pandas as pd
from itertools import repeat

sys.path.append("/hpc/local/CentOS7/dhl_ec/software/ctapp/")

from application.MongoParsenStore import reference_creator
from application.MongoParsenStore import liftover_creator
from application.processor import Processor

config_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'

funcs = {
        'liftover': liftover_creator.snp_liftover,
        '1000g': reference_creator.init_operations,
        'synonyms': liftover_creator.snp_synonyms,
        'GWAS_process': Processor.GWAS_preprocess
        }


class MapChunkers:

    def __init__(self, config, index):
        self.conf = config
        self.index = index
        self.doc = self.__load_config(config)
        self.data = self.doc.get('data_info')
        self.parameters = self.__param_check()
        self.__get_items()

    def __param_check(self):
        # list parameters may be a single item applicable to
        # all data files, in this case make an equal length iterable
        parameters = self.doc.get('parameters')
        if len(parameters) > 2:
            parameters = repeat(parameters, len(self.data))
        return parameters

    def __get_items(self):
        # get the wanted items from the config, based on given index
        self.params = self.parameters
        self.items = [self.data[ix] for ix in self.index]

    def __set_data_attr(self, doc):
        self.studyID = doc.get('studyID')
        self.path = doc.get('meta_data').get('path')
        doc = doc.get('properties')
        self.sep = doc.get('separator')
        self.skip = doc.get('skip')
        self.names = doc.get('headers')
        self.columns = doc.get('header_indices')
        if isinstance(self.columns, list):
            self.header_ix = self.skip
        else:
            self.header_ix = None

    def __set_param_attr(self, param):
        func = param.get('operation')
        self.func = funcs.get(func)
        self.mongo = param.get('mongo')
        self.chunksize = param.get('chunksize')
        self.subchunk = param.get("subchunk")
        if param.get("index_col"):
            self.index_col = param.get('index_col')
        else:
            self.index_col = None

    @staticmethod
    def __load_config(path):
        try:
            with open(path) as yaml_file:
                meta_yaml = yaml.load(yaml_file)
            return meta_yaml
        except yaml.YAMLError as exc:
            raise Exception("Not a valid YAML: {}".format_map(exc))
        except FileNotFoundError:
            raise Exception("No such file: {}".format(path))

    @staticmethod
    def __get_mongo_client():
        from pymongo import MongoClient

        def get_host_info():
            with open(config_path, 'r') as mon_conf:
                lines = mon_conf.readlines()
            return lines[1].split('=')[1].strip()

        monip = get_host_info()
        client = MongoClient(host=monip, connect=False)
        return client

    @staticmethod
    def __get_taskid():
        try:
            taskid = os.environ["SGE_TASK_ID"]
            return taskid
        except AttributeError:
            return False

    def __chunk_it(self, path, chsize, fskip, head_ix, sep, cols, ixcol, cnames, sub):
        # 1-indexed to 0-indexed
        taskid = self.__get_taskid()
        Nch = int(taskid)-1
        skip = Nch * chsize
        numrows = chsize
        numrows = None if numrows == 0 else numrows
        cols = None if len(cols) == 0 else cols
        dataframe_iter = pd.read_csv(filepath_or_buffer=path, sep=sep, header=head_ix, skiprows=skip, usecols=cols,
                                     index_col=ixcol, names=cnames, compression='infer', nrows=numrows, chunksize=sub,
                                     skipinitialspace=True)
        # return object is a dataframe iterator with the chunks as iter objects, if chunksize is not None
        return dataframe_iter

    def single(self):
        """
        single data processor, no data is not chunked so operation
        is only performed once
        :param index: index of the study_doc that is used from the config
        :return:
        """
        f_args = list()
        for item, param in zip(self.items, self.params):
            self.__set_data_attr(item)
            self.__set_param_attr(param)
            df = pd.read_csv(filepath_or_buffer=self.path, sep=self.sep, header=self.header_ix, usecols=self.columns,
                             index_col=self.index_col, names=self.names, skiprows=self.skip, compression='infer')
            f_args.append(df)
        if self.func is Processor.GWAS_preprocess:
            f_args.append(self.conf)
            f_args.append(self.index)
        if self.mongo:
            client = self.__get_mongo_client()
            f_args.append(client)
        # call function once for this job
        self.func(*tuple(f_args))

    def parallel(self):
        """
        parallel chunk processor, the function is executed with specific taskID related chunks
        :param index: index of the study_doc that is used from the config
        :return:
        """
        f_args = list()
        for item, param in zip(self.items, self.params):
            self.__set_data_attr(item)
            self.__set_param_attr(param)
            chunk_iterator = self.__chunk_it(self.path, self.chunksize, self.skip, self.header_ix, self.sep,
                                             self.columns, self.index_col, self.names, self.subchunk)
            f_args.append(chunk_iterator)
        if self.func is Processor.GWAS_preprocess:
            f_args.append(self.conf)
            f_args.append(self.index)
        if self.mongo:
            client = self.__get_mongo_client()
            f_args.append(client)
        # call function, with or without mongo client
        self.func(*tuple(f_args))

    def subparallel(self):
        """
        parallel chunk processor with subchunks as attribute, with use of threading
        :param index:
        :return:
        """
        from concurrent.futures import ThreadPoolExecutor
        threader = ThreadPoolExecutor(10)
        f_args = list()
        for item, param in zip(self.items, self.params):
            self.__set_data_attr(item)
            self.__set_param_attr(param)
            chunk_iterator = self.__chunk_it(self.path, self.chunksize, self.skip, self.header_ix, self.sep,
                                             self.columns, self.index_col, self.names, self.subchunk)
            if isinstance(chunk_iterator, pd.DataFrame):
                chunk_iterator = repeat(chunk_iterator)
            f_args.append(chunk_iterator)
        if self.func is Processor.GWAS_preprocess:
            f_args.append(repeat(self.conf))
            f_args.append(repeat(self.index))
        if self.mongo:
            client = self.__get_mongo_client()
            f_args.append(repeat(client))
        # call the functions in a parallel asynchronous fashion, with 'client' as a constant (when added)
        list(threader.map(self.func, *tuple(f_args)))  # launch functions simultaneously for every chunk in iterator


if __name__ == '__main__':
    fire.Fire(MapChunkers)
