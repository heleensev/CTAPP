import re, sys, os, logging, copy
from pymongo import ASCENDING
import platform
import timebuddy

# some notes: REF or ALT may contain multiple AA's separated by ',' eg: A,T means: multi allelic
# input file to process: 1000 Genomes Project phase3 final variant call set with 88 million variants
# site: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_phase3_callset_20150220

# logfile to log errors and significant events during runtime

logpath = '/hpc/local/CentOS7/dhl_ec/software/ctapp/logs'
sys.path.append('/hpc/local/CentOS7/dhl_ec/software/ctapp/')

path = os.path.join(logpath, 'reference_creator.log')
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

handler = logging.FileHandler(path)
handler.setLevel(logging.INFO)

formatter = logging.Formatter("%(threadName)s:%(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.info('starting reference_creator at: {}'.format(timebuddy.time()))

"""
snp doc example
{
    snp_id : rs3299203,
    chr : 1
    loc : 38922389
    ref : A
    qual : 100
    var : [
            {
                ALT: C
                AF: 0.23
                EUR_AF: 0.33
                AFR_AF: 0.55
            },
            {
                ALT: T
                AF: 0.36
                EUF_AF: 0.42
                AFR_AF: 015
            }
          ]
    }
"""


def init_operations(chunk, client):
    pipeline = process_rows(chunk)

    collection = init_db(client, dbname='1000g', colname='phase3_3')
    # create index for database on "SNP"
    collection.create_index([("SNP", ASCENDING)])
    bulk_insert(collection, pipeline, 'phase3')


def init_db(client, dbname, colname):
    db = client[dbname]
    collection = db[colname]
    return collection


# write pipeline to database
def bulk_insert(collection, pipeline, name):
    try:
        logger.info('len pipe {}= {}'.format(name, len(pipeline)))
        if len(pipeline) > 1:
            collection.insert_many(pipeline)
    except Exception as e:
        logger.error(e)


# process the values in the columns of the input chunk
def process_rows(chunk, snp_id=0):
    # new pipeline object for this chunk
    pipeline = list()
    try:
        iterable = chunk.itertuples(index=True)
        for i, chrom, loc, snp_id, ref, alt, qual, _, info in iterable:
            try:
                chrom, loc, snp_id, alt, alt_map, cqual = convert_type(chrom, loc, snp_id, ref, alt, qual)
                info_docs = process_info(ref, info, snp_id, alt, alt_map)

                for snp_doc in info_docs:
                    snp_doc['CHR'] = chrom
                    snp_doc['LOC'] = loc
                    snp_doc['REF'] = ref
                    snp_doc['QUAL'] = cqual
                    pipeline.append(snp_doc)
            except Exception as e:
                logger.error('error: {} in process rows'.format(e))
                logger.error('row:'.format(chrom, loc, snp_id, ref, alt, qual))
                logger.error('{} \n Error on line {}'.format(sys.exc_info()[1], sys.exc_info()[-1].tb_lineno))
                # skip this faulty row
    except:
        logger.error('{} \n Error on line {}'.format(sys.exc_info(), sys.exc_info()[-1].tb_lineno))
        logger.error('during process_columns at {}'.format(snp_id))
    return pipeline


def snp_div(alts, ref, snps):
    # snv: single nucleotide variant, dl: deletion, ins: insertion
    div_class = {'snv': [], 'dl': [], 'ins': []}
    divs = list()
    for a in alts:
        if len(a) == len(ref):
            div_class['snv'].append(a)
            divs.append('snv')
        elif len(a) < len(ref):
            div_class['dl'].append(a)
            divs.append('dl')
        else:
            div_class['ins'].append(a)
            divs.append('ins')
    dup = {}
    divs = [dup.setdefault(x, x) for x in divs if x not in dup]
    if len(divs) != len(snps):
        # if the number of different variants don't map to the total of snp ids, try to lump
        # insertions and deletions together, else map every alt allele to every snp id
        div = copy.deepcopy(div_class)
        indels = {'indel': [i for k in ('ins', 'dl') for i in div_class.get(k) if k]}
        if div_class.get('snv'):
            indels.update({'snv': div_class.get('snv')})
        div_class = indels
        divs = ['snv', 'indel'] if div.get('snv') else ['indel']
        if len(divs) != len(snps):
            # if no logical mapping can be done, map all the alleles to all the ids
            return {s: [i for a in div_class.values() for i in a] for s in snps}
    # make a map of the rs ids with the matching alt alleles
    div_to_snp = {div: snp for div, snp in zip(divs, snps)}
    snp_to_alt = {div_to_snp.get(s): v for s, v in div_class.items()}
    return snp_to_alt


def convert_type(chrom, loc, snp_id, ref, alt, qual):
    var_types = ['rs', 'esv', 'nsv', 'ss']
    try:
        # whole file= mixed datatypes, so condition does not hold
        chrom = str(chrom)
        if chrom.isnumeric():
            chrom = int(chrom)
        loc = int(loc)
        qual = int(qual)
        snps = snp_id.split(';')
        # characterizes dbvar snps
        if '<' not in alt:
            alts = alt.split(',')
        else:
            alts = [','.join(a.strip('<>').split(':')) for a in alt.split(',')]

        # if there is only one rs id and more than 1 alt allele, all are mapped to this id
        if len(snps) == 1 and len(alts) > 1:
            alt_map = {snps[0]: alts}
        # if not every alt allele can be mapped to a unique rs id, classify variant type
        elif len(snps) != len(alts):
            alt_map = snp_div(alts, ref, snps)
        # else every alt allele can be mapped to a unique rs id
        else:
            alt_map = {snp: [a] for a, snp in zip(alts, snps)}
        return chrom, loc, snps, alts, alt_map, qual
    except Exception as e:
        logger.error('error: {} in convert type'.format(sys.exc_info()))
        logger.error('snp: {}'.format(snp_id))
        logger.error('{} \n Error on line {}'.format(sys.exc_info()[1], sys.exc_info()[-1].tb_lineno))


def process_info(ref, info, snp_id, alts, alt_map):

    docs = list()
    info_dict = dict()
    aa_list = dict()
    flags = list()
    var_list = [{} for _ in alts]
    [var_list[x].update({'ALT': y}) for x, y in enumerate(alts)]
    info_dict['VAR'] = []
    aa_keys = ['ACA', 'RFA', 'ALA', 'IDT']
    # identifiers that may be present in the row, segmented in datatypes and formats
    info_units = [['CIEND', 'CIPOS', 'END', 'MEND', 'MLEN', 'MSTART', 'SVLEN', 'AC',
                   'NS', 'AN'], ['IMPRECISE', 'EX_TARGET', 'MULTI_ALLELIC'], ['AA'],
                  ['AF', 'EAS_AF', 'AFR_AF', 'AMR_AF', 'EUR_AF', 'SAS_AF'],
                  ['CS', 'MC', 'MEINFO', 'SVTYPE', 'OLD_VARIANT']]
    # corresponding patterns for the identifiers, with lookahead & lookbehind to capture the wanted characters
    expressions = ['(?<={}=)\d+(?=;|$)', '{}', '(?<={}=)[\w\|\.\(\)-]*(?=;|$)',
                   '(?<={}=)[(\d)\.,]+(?=;|$)', '(?<={}=)[-_(\w),]+(?=;|$)']
    try:
        # simultaneously loop over ids in info_units and expression in expressions
        for loop, (unit, expression) in enumerate(zip(info_units, expressions)):
            for type_id in unit:
                pattern = re.compile(r'{}'.format(expression.format(type_id)), re.M | re.I)
                search = pattern.search(info)
                if search:
                    match = search.group()
                    # depending on loop number, process match
                    if loop == 0:
                        val = int(match)
                    elif loop == 1:
                        flags.append(match)
                        continue
                    elif loop == 2:
                        aas = match.split('|')
                        aas = [x if not re.match(r'(\.)|(unknown.*)', x) else '' for x in aas]
                        [aa_list.update({k: v}) for k, v in zip(aa_keys, aas)]
                        val = aa_list
                    elif loop == 3:
                        var_frqs = [float(x) for x in match.split(',')]
                        for x, y in enumerate(var_frqs):
                            var_list[x].update({type_id: y})
                        continue
                    else:
                        val = str(match)
                    # append the info as key, val to info_list
                    info_dict[type_id] = val
        info_dict['flags'] = flags
        # use alt_map to map the alt alleles to the right SNP identifier and
        # append the specific info of the variant to a separate dictionary
        for snp in alt_map:
            doc = copy.deepcopy(info_dict)
            for alt in alt_map.get(snp):
                i = alts.index(alt)
                # check the variant type of this allele
                if 'SVTYPE' in doc:
                    VT = 'SV'
                else:
                    VT = 'SNP' if len(alt) == len(ref) else 'INDEL'
                var_list[i].update({'VT': VT})
                doc['VAR'].append(var_list[i])
                doc['SNP'] = snp
            docs.append(doc)

        return docs
    except Exception as e:
        logger.error('error:{} in process_info'.format(e))
        logger.error('row: {} {} {}\ndoc: {}'.format(info, snp_id, alts, docs))
        logger.error('{} \n Error on line {}'.format(sys.exc_info()[1], sys.exc_info()[-1].tb_lineno))