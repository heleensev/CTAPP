import pandas as pd
import numpy as np
import logging
import math
import re
import sys

sys.path.extend('/hpc/local/CentOS7/dhl_ec/software/ctapp')

from application.GWASParser import snp_lookup

logging.basicConfig(filename='../logs/valchecker.log', filemode='a', level=logging.INFO)
logger = logging.getLogger()

col_pttrs = {'SNP': '((chr)?\d{1,2}[:;_-](\d)+)|(^(rs[ _-]?)?\d+)',
             'CHR': '(\d){1,2}|[XYMT]{1,2}',
             'BP': '\d{1,9}',
             'A1': '[ATCGRDI]+',
             'A2': '[ATCGRDI]+',
             'FRQ': '[\d.]+(?:e-?\d+)?',
             'BETA': '-?[\d.]+(?:e-?\d+)?',
             'OR': '-?[\d.]+(?:e-?\d+)?',
             'SE': '-?[\d.]+(?:e-?\d+)?',
             'P': '-?[\d.]+(?:e-?\d+)?',
             'N': '(\d)+',
             'INFO': '-?[\d.]+(?:e-?\d+)?'}


def check_correct(data, study):

    def snp_check(v, regex):
        val = str(v)
        match = regex.match(val)
        if match:
            # check if rs prefix is present, if not, return concatenated value
            if match.group(4) and match.group(5):
                return v
            elif match.group(4):
                return 'rs' + match.group()
            elif match.group(1):
                if not match.group(2):
                    val = 'chr{}:{}'.format(val.split(':')[0], val.split(':')[1])
                # call snp_lookup
                rs = snp_lookup.liftover_loc(study, val)
                if rs:
                    return rs
                else:
                    return val
        print(v)
        return np.NaN

    def chr_check(v, regex):
        val = str(v)
        if regex.match(val):
            # NCBI notation for X, Y, XY and MT
            chr_map = {'X': 23, 'Y': 24, 'XY': 25, 'MT': 26}
            if val in chr_map:
                return chr_map.get(val)
            # check if autosomal chromosome number is no higher than 26
            elif val.isnumeric() and int(val) < 27:
                return v
        return np.NaN

    def bp_check(v, regex):
        val = int(v)
        if regex.match(str(val)):
            bp_human_genome = 300000000
            if val < bp_human_genome:
                return v
        return np.NaN

    def all_check(v, regex):
        val = str(v)
        if regex.match(v, re.I):
            # old PLINK notation for alleles
            plink_map = {'1': 'A', '2': 'C', '3': 'G', '4': 'T'}
            # multiple allele check?
            if val in plink_map:
                return plink_map.get(val)
            elif val.islower():
                return val.upper()
            return v
        return np.NaN

    def frq_check(v, regex):
        val = str(v)
        # convert scientific format
        if regex.match(val):
            return v
        return np.NaN

    def eff_check(v, regex):
        val = str(v)
        if regex.match(val):
            if study.effect_type == 'OR':
                # return log odds (beta) if value is odds ratio
                # if value is negative, return ln of positive value
                try:
                    math.log(v)
                    return math.log(v)
                except ValueError:
                    return math.log(v*-1)
            else:
                return v
        return np.NaN

    def se_check(v, regex):
        val = str(v)
        if regex.match(val):
            return v
        return np.NaN

    def p_check(v, regex):
        val = str(v)
        if regex.match(val):
            return v
        return np.NaN

    def n_check(v, regex):
        val = str(v)
        if regex.match(val):
            if v < 0.9:
                return np.NaN
            return v
        return np.NaN

    def info_check(v, regex):
        val = str(v)
        if regex.match(val):
            return v
        return np.NaN

    def pattern(head):
        pattern = col_pttrs.get(head)
        return re.compile(r'{}'.format(pattern), re.I)

    columns = ['SNP', 'CHR', 'BP', 'A1', 'A2', 'FRQ', 'SE', 'P', 'BETA', 'OR', 'N', 'INFO']

    col_types = [snp_check, chr_check, bp_check, all_check, all_check, frq_check,
                 se_check, p_check, eff_check, eff_check, n_check, info_check]

    data_headers = list(data.columns.values)

    for col, ctype in zip(columns, col_types):
        if col in data_headers:
            data[col] = data[col].apply(lambda x: ctype(x, pattern(col)))
    if 'OR' in data_headers:
        data_headers[data_headers.index('OR')] = 'BETA'
        data.columns = data_headers

    # obtain SNP with BED format as ID instead of rs id
    rs_mask = data['SNP'].str.contains('chr')
    rs_mask.fillna(False, inplace=True)
    no_rs = data[rs_mask]
    rs = data[~rs_mask]

    # get SNPs of all dropped columns
    all_dropped = data['SNP'][data.isnull().any(axis=1)]
    # write all_dropped and nulls to logfile
    logger.info('all dropped SNPS because of null values in {0}\n{1}'
                .format(study.studyID, all_dropped))

    # drop rows containing NaN's
    data.dropna(inplace=True)

    return rs, no_rs

def test_val_checker():
    class study:
        def __init__(self, eff):
            self.effect = eff

    df = pd.read_csv('/home/sevvy/Scripts/deploy/data/test/test_s.csv', header=0, sep='\t')
    stu = study('OR')
    check_correct(df, stu)

{ "VAR" : [ { "ALT" : "CTTTTTTTT", "VT" : "INDEL" }, { "ALT" : "CTTTTTTTTT", "VT" : "INDEL" },
                                                     { "ALT" : "CTTTTTTTTTT", "VT" : "INDEL" },
                                                     { "ALT" : "CTTTTTTTTTTC", "VT" : "INDEL" },
                                                     { "ALT" : "CTTTTTTTTTTT", "VT" : "INDEL" } ], "NS" : 2504, "AN" : 5008, "flags" :
    [ "MULTI_ALLELIC" ], "SNP" : "rs5847627", "CHR" : 3, "LOC" : 30797316, "REF" : "C", "QUAL" : 100 }




