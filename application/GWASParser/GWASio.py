import logging
import re
import os
from glob import glob
import pandas as pd

filename = str()
logger = logging.getLogger(__name__)

magma_datadir = '/hpc/local/CentOS7/dhl_ec/software/ctapp/data/magma'


def init_reader(study):
    logger.info("Entering read_GWAS")
    path = study.path
    sep = check_sep(path)
    # update studyID doc with separator
    study.sep = sep

    return study


def open_file(file):
    global filename
    try:
        # read first line of file, if line is longer than 1000, probably not a
        # header line or line was not proper split because of lacking new lines
        with open(file) as GWAS_file:
            header = GWAS_file.readline(1000)
            if header:
                return header
            else:
                if '\n' not in header:
                    raise Exception("file does not contain newlines, "
                                    "or has extremely long lines")
    except Exception as e:
        logger.error("Error occurred during separator check: {}".format(e))


def check_sep(path):
    try:
        # call open_file to check for a valid separator, and a valid header line
        header_line = open_file(path)
        logger.info("entering check sep\n")
        b = {}
        # find all non-newlines in the header line
        separators = re.findall(r"\W", header_line)
        # for all found non-alphanumeric & non-newlines, write to a dictionary
        # and count entries by adding up previous identical entries
        for sep in separators:
            b[sep] = b.get(sep, 0) + 1
        # list for
        cnts = list()
        seps = list()
        for sep in b.keys():
            seps.append(sep)
            cnt = b.get(sep)
            cnts.append(int(cnt))
        max_cnt = max(cnts)
        if max_cnt > 5:
            sep = seps[cnts.index(max_cnt)]
            logger.info("valid separator found: " + sep + "\n")
            return sep
        else:
            mssg = "no valid separator found"
            print(mssg)
            raise Exception(mssg)
    except Exception as e:
        logger.error("Error occured during separator check: {}".format(e))


def GWAS_out(data, ID):

    """
    Final step of GWAS preprocessor, prepare GWAS csv's for MAGMA processing
    by writing the checked and corrected data to a new csv file. The columns
    are sorted to be conform to the desired MAGMA format: SNPid, CHR, BP, P, (N).
    :param data:
    :param ID:
    :return:
    """
    # import threading for identification of current thread
    import threading
    taskid = os.environ["SGE_TASK_ID"]

    rs, no_rs = data

    for data, prefix in zip([rs, no_rs], ['rs', 'no_rs']):
        if not data.empty:
            cols = list(data.columns)

            # list unwanted columns
            to_drop = [cols.index(c) for c in cols if c not in ['SNP', 'CHR', 'BP', 'P', 'N']]
            # drop unwanted columns
            data.drop(data.columns[[to_drop]], axis=1, inplace=True)

            # current columns in data
            cols = list(data.columns)

            # reorder columns, put N column at the end if present
            ordered_cols = ['SNP', 'CHR', 'BP', 'P']
            if 'N' in cols:
                ordered_cols += ['N']

            data = data.reindex(columns=ordered_cols, copy=False)
            data.columns = ordered_cols

            # get thread id to discern chunk from other processes in the filename
            thr = re.split('[-_,]', str(threading.current_thread()))[2]

            filename = '{}_{}_{}_{}.csv'.format(prefix, ID, taskid, thr)
            filepath = os.path.join(magma_datadir, filename)

            data.to_csv(filepath, sep='\t', index=None, header=True)


def concat(ID):
    """
    concatenator for all the processed chunks, as preparation for MAGMA
    gene annotation and gene analysis
    :return:
    """
    from glob import glob

    for prefix in ['rs', 'no_rs']:
        filenames = sorted(glob(os.path.join(magma_datadir, '{}_{}_*_*.csv'.format(prefix, ID))))
        if filenames:
            complete_df = pd.concat((pd.read_csv(f, header=0, sep='\t', dtype={'CHR': int, 'BP': int}) for f in filenames))
            filename = os.path.join(magma_datadir, '{}_{}_complete.csv'.format(prefix, ID))
            complete_df.to_csv(filename, header=True, sep='\t', index=False)


def merge_chunks(studyID):

    filenames = sorted(glob(os.path.join(magma_datadir, studyID + '*.csv')))

    complete_df = pd.concat((pd.read_csv(f, header=0, sep='\t') for f in filenames))

    filename = os.path.join(magma_datadir, '{}_complete.csv'.format(studyID))
    complete_df.to_csv(filename, header=True, sep='\t', index=False)

