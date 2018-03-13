from timebuddy import time
from pymongo import errors
import numpy as np
import logging, sys

sys.path.append('/hpc/local/CentOS7/dhl_ec/software/ctapp')
from application.GWASParser import snp_lookup

frq_switch_sum = 0
no_match_sum = 0

logging.basicConfig(filename='../logs/refchecker.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger(__name__)


def init_db(client, dbname, colname):
    db = client[dbname]
    collection = db[colname]
    return collection


def init_ref_check(data, study, client):
    rs, no_rs = data
    pop = str(study.ethnicity).upper() + '_AF'
    ref_doc = fetch_ref_alleles(rs, client, pop)

    cols = ['SNP', 'BP', 'EA', 'OA', 'FRQ', 'BETA']
    # coerce float64 so datatype is compatible with numpy NaN
    num_cols = ['BP', 'FRQ', 'BETA']
    rs[num_cols] = rs[num_cols].astype('float64')

    rs[cols] = rs[cols].apply(reference_check, axis=1, args=(ref_doc, pop))
    # if dataframe contains SNPs with no rs id, perform separate reference check
    if not no_rs.empty:
        no_rs[cols] = no_rs[cols].apply(ref_check_no_rs, axis=1, args=study.build)

    return rs, no_rs
    # https://stackoverflow.com/questions/40353519/how-to-apply-custom-function-to-pandas-data-frame-for-each-row
    # swap around by returning different order


def ref_check_no_rs(x, build):
    """
    reference check for the SNPs with no 'rs' identifier, alleles and corresponding values must
    be mapped to a standard format to make meta anlysis possible
    :param x:
    :return:
    """
    SNP, BP, EA, OA, FRQ, BETA = x['SNP'], x['BP'], x['EA'], x['OA'], x['FRQ'], x['BETA']

    # replace the old bp positions with hg19 if it is not original build
    if build != 'hg19':
        x['BP'] = snp_lookup.liftover_to_19(x['BP'], build)

    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # arbitrary reference strand sequence nucleotide start
    if EA.startswith('T') or EA.startswith('C'):
        x['EA'] = ''.join([complements.get(nuc) for nuc in EA])
        x['OA'] = ''.join([complements.get(nuc) for nuc in OA])
        EA = x['EA']
        OA = x['OA']
    # sort the alleles on length (ie indels become the 'alt' allele) and alphabetic order
    alleles = [EA, OA]
    if len(EA) > len(OA) or not alleles == sorted(alleles):
        x['EA'], x['OA'] = x['OA'], x['EA']

        x['FRQ'] = 1-FRQ
        x['BETA'] = BETA*-1

    return x


def fetch_ref_alleles(data, client, pop):
    def bulk_find():
        try:
            results = collection.find({"SNP": {"$in": id_array}}, {"_id": 0, "LOC": 1, "REF": 1, "VAR": 1, "SNP": 1})
            return results
        except errors.ConnectionFailure:
            print('ehhh')

    collection = init_db(client, '1000g', 'phase3')
    # turn the column containing the rs numbers into a list
    id_array = data['SNP'].tolist()
    # find corresponding rs ids in reference (1000 genome) database
    docs = list(bulk_find())
    # refs = [[d['REF'], d['LOC']] for d in docs]
    # get all the variant info (from the results) per snp into one list
    vars = [d['VAR'] for d in docs]
    refs = [d['REF'] for d in docs]
    locs = [d['LOC'] for d in docs]
    # get all the rs ids (from the results) into one list
    snps = [d['SNP'] for d in docs]
    # put the corresponding rs id and the wanted variant info together in one dictionary
    var_doc = [[{k: v[k] for k in (pop, 'ALT')}] for snp, var in zip(snps, vars) for v in var]
    doc = {snp: {'LOC': loc, 'REF': ref, 'VAR': var} for loc, snp, var, ref in zip(locs, snps, var_doc, refs)}
    ref_doc = doc
    return ref_doc


def reference_check(x, ref_doc, pop):
    # check the SNPs in the GWAS set with the reference (1000 genome call set)
    # desired format: A1: ref allele, EF: alt allele, FRQ: alt allele frequency
    names = ['SNP', 'BP', 'EA', 'OA', 'FRQ', 'BETA']
    SNP, BP, EA, OA, FRQ, BETA = x['SNP'], x['BP'], x['EA'], x['OA'], x['FRQ'], x['BETA']

    # in orientation check the allele corresponding to the reference allele is found
    # bool for (not) removing row at the end of this function block
    try:

        def orientation_check(firstcheck=True, swapfirst=True):
            # allele 1 and 2, may be changed throughout this function block,
            # if orientation differs from reference
            # complementary alleles to look up during orientation check
            complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            global EA, OA, swap
            if swapfirst:
                # allele_check if the ref and alt allele corresponds with assumed position
                if ref_all == x['OA'] and alt_all == x['EA']:
                    # wanted orientation, so alleles don't need swapping
                    swap = False
                    return True
                # allele_check if ref and alt allele are not oriented to reference
                elif ref_all == x['EA'] and alt_all == x['OA']:
                    # swap the two variables around
                    x['EA'], x['OA'] = x['OA'], x['EA']
                    swap = True
                    return True
            # if this check was already performed, function will return False
            elif firstcheck:
                # # if alleles are complementary to the reference, these are altered
                x['EA'] = ''.join([complements.get(a) for a in EA])
                x['OA'] = ''.join([complements.get(a) for a in OA])
                if orientation_check(firstcheck=False):
                    logger.info('complementary alleles for {}'.format(SNP))
                    return True

        def frequency_check():
            global frq_switch_sum, swap
            frq_range = 0.15
            eur_frq = float(var.get(pop))

            def check_range(inverse=False):
                FRQ = x['FRQ']
                FRQ = float('%f' % FRQ)
                if inverse:
                    FRQ = 1 - FRQ
                under_range = FRQ - frq_range
                upper_range = FRQ + frq_range
                # if eur_freq is smaller than under_range and higher than upper_range
                if under_range < eur_frq < upper_range:
                    x['FRQ'] = FRQ
                    x['BETA'] = (x['BETA'] * -1)
                    return True

            # check most frequently occurring condition first
            if check_range(inverse=swap):
                # if the alleles were swapped, frequency and odds ratio need to be transformed
                if swap:
                    frq_switch_sum += 1
                return True
            elif check_range(inverse=not swap):
                # if the swapped alleles do have matching frequency to te reference
                # may be due to A/T C/G SNPs, flip them back and perform orientation check again
                x['EA'], x['OA'] = x['OA'], x['EA']
                if orientation_check(swapfirst=False):
                    frq_switch_sum += 1
                    return True

        def indel_check():
            global swap
            ref_variant = var.get('VT')
            if ref_variant == 'INDEL':
                if x['OA'] in ['RDI']:
                    swap = True
                x['OA'] = ref_all
                x['EA'] = var.get('ALT')
                return True

        # calling of the above functions

        # search database for matching SNP identifier
        lookup = ref_doc.get(SNP)
        if not lookup:
            # return NaN for all the columns
            return np.NaN
        # for multi_allelic SNPs, multiple alleles for one rs ID, so all are checked
        ref_all = lookup.get('REF')
        vars = lookup.get('VAR')
        for i, var in enumerate(vars):
            swap = False
            alt_all = var.get('ALT')
            if OA in ['RDI'] or EA in ['RDI']:
                if not indel_check():
                    continue
            else:
                if not orientation_check():
                    continue
            if frequency_check():
                # position in reference dataset (GRCh37) may differ from GWAS dataset
                ref_loc = lookup.get('LOC')
                # replace base pair position if different from reference
                if int(x['BP']) != ref_loc:
                    x['BP'] = ref_loc
                return x
            else:
                continue
                # or np.NaN?
        # return NaN for all values in row if no match has been found
        for i in names:
            x[i] = np.NaN
        return x

    except Exception as message:
        logger.info('{} for SNP: {} ({},{})'.format(message, SNP, EA, OA))
    except:
        logger.error(sys.exc_info())
        logger.error('during liftover at row {}'.format(SNP))


"""
check order unpacking of values from iterable after reading file
rs ID may not be present in the file, so no liftover, but alleles and frequency may be checked
turned cursor into list at fetching of documents
"""
