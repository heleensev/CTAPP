# Module for automatic check on datatypes in the individual columns
import re
import sys
import logging
import traceback
import pandas as pd

# sys.path.append('/home/sevvy/Scripts/deploy/application')
# from application.GWASParser import GWASio
#
from . import GWASio

logger = logging.getLogger(__name__)

col_types = [['SNP', '(snp)|(marker[ -_]?(name)?)|(rs[ _-]?(id))', '((rs[ _-]?)?\d+)|'
              '((chr)?\d{1,2}\:(\d)+)'],
             ['CHR', '(chr(omosome)?)', '[1-26]|[XY]'],
             ['BP', '((pos\Z)|(loc(ation)?)|(bp\Z)|(hg(\d){2})|(grch(\d){2}))', '\d+'],
             ['Allele', '((allele)|(A)[-_ ]?[12]?)|major|minor|(non[ _-]?)?effect|other|risk',
              '[ACTG]+|[RDI]?|[1234]?'],
             ['FRQ', 'fr(e)?q|maf|eaf', '(\d)*((\.)(\d)*(E)?(-)?(\d)*)?'],
             ['EFF', '(beta)|(effect)|(OR)|odds[ _-]?ratio', '(-)?(\d)*((\.)(\d)*(E)?(-)?(\d)*)?'],
             ['SE', '(se)|(std(\w)*)|((standard( -_)?)?error)', '(\d)*((\.)(\d)*(E)?(-)?(\d)*)?'],
             ['N', '(^N([ _-])?)?((studies)|(case(s)?)|(control(s)*)|$)', '[0-1000000]+'],
             ['P', 'p([ _-])?\.?(val)?(ue)?', '(-)?(\d)*((\.)(\d)*(E)?(-)?(\d)*)?'],
             ['INFO', '(info)|(imputation)|(variance)', '(\d)*((\.)(\d)*(E)?(-)?(\d)*)?']]


def init_classifier(study):
    this_study = GWASio.init_reader(study)
    header_IDer(this_study)
    post_checks(this_study)

    return this_study


def header_IDer(this_study):
    path = this_study.path
    sep = this_study.sep

    def allele_check():
        new_header = dup_vals_check('Allele', ['A1', 'A2'])
        if new_header:
            return new_header

    # check and correction for values that are allowed to be duplicate in the headers (i.e: allele)
    def dup_vals_check(val, dup_vals):
        # if the current column name from col_patterns in the loop matches Allele
        if col[0] == val:
            if dup_vals[0] not in new_headers:
                new_head = dup_vals[0]
            elif dup_vals[1] not in new_headers:
                new_head = dup_vals[1]
            else:
                raise Exception('duplicate header found: {}'.format(header))
            return new_head

    def duplicate_check():
        # if header is already in the headers list (for header in headers list except current header)
        # headers minus current evaluated value
        headers_min_cur = [x for n, x in enumerate(new_headers) if n != i]
        if col[0] in headers_min_cur:
            this_study.err_mssg = 'duplicate header: {}'.format(new_header)
        return True

    # transform the csv to a DataFrame, to determine the headers
    # get a chunk of the file to perform the check (for memory efficiency)
    df = pd.read_csv(filepath_or_buffer=path, sep=sep, nrows=25, header=0)
    headers = df.columns.values.tolist()
    # add original headers as attribute to the study object
    this_study.org_headers = headers

    dispose = list()
    new_headers = list()
    header_idx = list()

    # loop over all columns to determine the information type
    for i, header in enumerate(headers):
        try:
            df_chunk = df[header]
            for c, col in enumerate(col_types):
                if col_check(df_chunk, header, this_study, col[1], col[2]):
                    new_header = col[0]
                    # call duplicate check for duplicate headers
                    if new_header == 'Allele':
                        new_header = allele_check()
                        new_headers.append(new_header)
                        header_idx.append(i)
                        break
                    elif duplicate_check():
                        # replace old header in list with new header
                        new_headers.append(new_header)
                        header_idx.append(i)
                        # break loop when column is identified
                        break
            else:
                # if no valid header match found for the column, dispose column
                dispose.append(header)
        except Exception as e:
            if e == 'duplicate':
                this_study.success = False
                continue
            print(traceback.format_exc())
            logger.error(traceback.format_exc())

    this_study.disposed = ' '.join(dispose)
    # set new attributes to the study object
    this_study.headers = new_headers
    this_study.head_idx = header_idx


def col_check(df, header, study, re_head, re_row):
    head, row = False, False
    hd_pattern = re.compile(r'{}'.format(re_head), re.I)
    row_pattern = re.compile(r'{}'.format(re_row), re.I)

    match = [row_pattern.match(str(val)) for i, val in df.iteritems()]
    if not None in match:
       row = True
    # if the column header matches the header pattern, the header is confirmed
    if hd_pattern.match(header):
        # header heterogeneity test values often contain 'het', this column is disposed of
        if not re.match(r'het', header, re.I):
            head = True
    # # if the type is confirmed, but not the header, header may actually  be a row value
    # # implies header is missing from the input file, call identical_increment to count occurrences
    # if row and row_pattern.match(header):
    #     study.err_mssg = 'Header missing for column, continuing to user_classify_columns'
    # if header or column matches with the pattern, header is confirmed
    if head and row:
        return True
    return False


def post_checks(study):

    required = ['SNP', 'CHR', 'BP', 'A1',
                'A2', 'FRQ', 'EFF', 'SE', 'P']
    # column with study sizes required if n_studies is none
    if not hasattr(study, 'n_studies'):
        required.append('N')

    headers = study.headers
    print(headers)
    header_idx = study.head_idx
    org_headers = study.org_headers
    effect_allele = study.effect_allele
    effect_type = study.effect_type

    def check_essential():
        reqs = set(required)
        heads = set(headers)
        missing = reqs.difference(heads)
        if missing:
            study.missing = list(missing)
            study.success = False

    def header_indices():
        allowed = ['SNP', 'CHR', 'BP', 'EA', 'OA', 'FRQ', 'EFF', 'SE', 'P', 'INFO', 'N']
        # get header indices mapped to the reference (required list)
        idx_to_ref = [allowed.index(head) for head in headers]
        study.head_to_ref = idx_to_ref

    def effect_allele_check():
        if effect_allele not in org_headers:
            study.err_mssg = 'ERROR: effect allele name not in original headers'

        ea_header_idx = org_headers.index(effect_allele)
        allele_idx = [headers.index('A1')] + [headers.index('A2')]
        real_idx = [header_idx[allele_idx[0]]] + [header_idx[allele_idx[1]]]
        if ea_header_idx not in real_idx:
            study.err_mssg = 'ERROR: effect allele column is not identified as allele column'
        else:
            headers[ea_header_idx] = 'EA'
            NFA_ix = allele_idx[0] if allele_idx[0] != ea_header_idx else allele_idx[1]
            headers[NFA_ix] = 'OA'
        if study.err_mssg:
            study.success = False

    def replace_effect_type():
        eff_idx = headers.index('EFF')
        headers[eff_idx] = effect_type

    check_essential()
    effect_allele_check()
    header_indices()
    replace_effect_type()

# class Study:
#     def __init__(self):
#         self.path = '/home/sevvy/Scripts/deploy/data/test/val_test.csv'
#         self.sep = '\t'
#         self.effect_allele = 'A1'
#         self.effect_type = 'BETA'
#         self.err_mssg = ''
#
# s = Study()
# init_classifier(s)
