import fire
import os
import logging
import timebuddy
import pandas as pd
from pymongo import UpdateOne, ASCENDING, errors

logging.basicConfig(filename='../logs/genedb.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger(__name__)

monconf_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'
magma_io = '/hpc/local/CentOS7/dhl_ec/software/ctapp/data/MAGMA'


def get_mongo_client():
    from pymongo import MongoClient

    def get_host_info():
        with open(monconf_path, 'r') as mon_conf:
            lines = mon_conf.readlines()
        return lines[1].split('=')[1].strip()

    monip = get_host_info()
    client = MongoClient(host=monip, connect=False)
    return client


def init_db(client, dbname, colname):
    db = client[dbname]
    col_names = list(db.collection_names())
    # if collection already exists, append suffix
    # if colname in col_names:
    #     colname = '{}_1'.format(colname)
    collection = db[colname]
    return collection


def bulk_exec(col, bulk_obj):
    try:
        col.bulk_write(bulk_obj, ordered=True)
    except Exception as e:
        pass


def gene_update(studyID, prefix, annal):
    """
    update the genes database for the GWAS results
    :param studyID: ID of GWAS study
    :param prefix: prefix of MAGMA out file to process
    :return:
    """

    client = get_mongo_client()
    colname = ''.join(timebuddy.date().split('-'))
    gene_col = init_db(client, 'MetaGWAS', 'genes{}'.format(colname))
    # create indexes if it not already exists
    gene_col.create_index([('GENE', ASCENDING)], unique=True)
    gene_col.create_index([('STUDY.StudyID', ASCENDING)])

    # get SNPS collection for updating GENE reference
    snps_col = init_db(client, 'MetaGWAS', 'snps{}'.format(colname))

    # get geneID collection, for annotation of uniprot ids
    geneID_col = init_db(client, 'IDmapping', 'hgnc_ids2')

    # parse MAGMA gene annotation output and write to database
    parse_MAGMA_annot(studyID, prefix, gene_col, snps_col, geneID_col)

    # if only the annotation was done, terminate here
    if prefix == 'no_rs' or not annal:
        return
    # if MAGMA gene analysis was done, continue to gene analysis parser
    parse_MAGMA_genes(studyID, gene_col)


def collect_ids(genes_list, col):
    try:
        id_annot = list(col.find({'entrez': {'$in': genes_list}}))

        annot_dict = dict()
        for gene_doc in id_annot:
            if gene_doc['entrez']:
                annot_dict[gene_doc['entrez']] = gene_doc
        return annot_dict

    except errors.BulkWriteError:
        logger.info("")
    except Exception as e:
        logger.error("")


def parse_MAGMA_annot(ID, prefix, gene_col, snp_col, geneID_col):
    annot_file = os.path.join(magma_io, '{}_{}.genes.annot'.format(prefix, ID))

    bulk = list()
    # nps bulk update list
    snp_bulk= list()
    # read all genes from the MAGMA gene annotation file
    all_genes = pd.read_csv(annot_file, skiprows=2, usecols=[0], sep='\t', header=None, squeeze=True)
    genes_list = list(all_genes)
    # get gene annotation (mapping to multiple gene ids) from geneID db, function returns a dictionary
    genes_annot = collect_ids(genes_list, geneID_col)

    with open(annot_file, 'r') as MAGMAannot:
        # read first two headers, ie remove first two lines from generator
        for _ in range(2):
            next(MAGMAannot)
        # get annotation info per line
        for line in MAGMAannot:
            line = line.split('\t')
            gene = line[0]
            loc = line[1].split(':')
            chr, start, stop = (loc[0], loc[1], loc[2])
            snps = line[2:]
            snps[-1] = snps[-1].strip('\n')

            # get additional info from geneID db
            gene_doc = genes_annot.get(gene)
            if gene_doc:
                annot = {gene_doc['uniprot'], gene_doc['symbol'], gene_doc['name'], gene_doc['ensembl']}

            # query filter, creates document if it does not exist already
            gene_filter = {'GENE': gene, 'CHR': chr, 'START': start, 'STOP': stop}
            # if present, add gene annotation
            if annot:
                gene_filter.update(annot)
            # modify the update query so that it contains the study ID and the push operator
            gene_update = {'$push': {'STUDY': {'StudyID': ID, 'SNPS': snps}}}

            # for the snp document
            # query filter selects snps associated with the gene
            snp_filter = {'SNP': {'$in': snps}}
            # query update, adds a reference to the gene document
            snp_update = {'$set': {'GENE': gene}}

            # add the query update to the bulk pipelineresults
            bulk.append(UpdateOne(gene_filter, gene_update, upsert=True))
            # add the update query to the snp bulk pipeline
            snp_bulk.append(UpdateOne(snp_filter, snp_update))

    bulk_exec(gene_col, bulk)
    bulk_exec(snp_col, bulk)


def parse_MAGMA_genes(ID, col):

    genes_file = os.path.join(magma_io, '{}.genes.out'.format(ID))
    genesdf = pd.read_csv(genes_file, delim_whitespace=True, header=0)

    # gene ids function as filter
    geneids = genesdf['GENE'].astype(str)
    # gene parameters that are study specific
    studyparam = ['NSNPS', 'NPARAM', 'N', 'ZSTAT', 'P', 'RSQ']
    studyspecs = genesdf[[studyparam]]
    # convert appropriate fields to floats
    studyspecs[['NPRAM', 'ZSTAT', 'P', 'RSQ']] = studyspecs[['NPRAM', 'ZSTAT' 'P']].astype(float)
    # convert appropriate fields to ints
    studyspecs[['NSNPS', 'N']] = studyspecs[['NSNPS', 'N']].astype(int)
    # turn dataframes into dictionaries to process in db insertion
    filters = geneids.to_dict(orient='records')
    updates = studyspecs.to_dict(orient='records')

    # genes bulk doc update list
    gene_bulk = list()

    for query_filter, query_update in zip(filters, updates):
        # modify query filter part of the query so that the update pushes to the right subdoc
        query_filter = query_filter.update({'STUDY': {'$elemMatch': {'StudyID': ID}}})
        # modify the update query so that it contains the push operator and updates on the filter position at $elemMatch
        query_update = {'$set': {'STUDY.$': query_update}}
        # add the update query to the gene bulk pipeline
        gene_bulk.append(UpdateOne(query_filter, query_update, upsert=True))

    bulk_exec(col, gene_bulk)

    "new function"

    genes_file = os.path.join(magma_io, '{}.genes.out'.format(ID))
    genesdf = pd.read_csv(genes_file, delim_whitespace=True, header=0)



    """
    db.col.update({"SNP": 12565280, "CHR": 1, "BP": 711153, "REF": "A"},
    {$push: {VAR: {"ALT": "C", "FRQ": 0.0662, "OR": -0.914, "SE": 0.38923}}}, {upsert:true})
    """

    """
    reqs = UpdateOne(filter={"SNP": 12565280, "CHR": 1, "BP": 711153, "REF": "A"},
    update={'$push': {'VAR': {"ALT": "CG", "FRQ": 0.4378, "OR": -0.380, "SE": 0.4785}}}, upsert=True)

    db.col.update({"GENE": 2398230, "CHR": 1, "START": 27387283, "STOP": 27388000},
                    {$push: {STUDY: {"NSNPS": 12, "NPARAM": 3, "N": 60, "ZSTAT": 0.45, "P": 0.003}}}, {upsert=true})

    # second, expand query to match the previously upserted alt allele
    # and add the study specific metadata of the SNP
    q_filter['VAR'] = {'$elemMatch': {'ALT': a2}}
    study_update = {'$push': {"VAR.$.STUDY": {snp_info}}}
                    """

def meta_update(chrno):
    """
    update gene database with MAGMA meta analysis results
    :param chrno: chromosome number, identifier for MAGMA output chunks
    :return:
    """
    client = get_mongo_client()
    date = ''.join(timebuddy.date().split('-'))
    col = init_db(client, 'GWASgenes', 'genes{}'.format(date))
    parse_MAGMA_genes(col, chrno)


def parse_MAGMA_meta(col, chrno):
    genes_file = os.path.join(magma_io, 'meta_{}.genes.out'.format(chrno))

    genesdf = pd.read_csv(genes_file, delim_whitespace=True, header=0)

    # gene ids function as filter
    geneids = genesdf['GENE'].astype(str)

    # gene parameters that are specific for the meta analysis
    studyparam = ['NSNPS', 'NPARAM', 'N', 'ZSTAT', 'P', 'DATASETS']
    studyspecs = genesdf[[studyparam]]
    # convert appropriate fields to ints
    studyspecs[['NSNPS', 'N', 'DATASETS']] = studyspecs[['NSNPS', 'N', 'DATASETS']].astype(int)
    # convert appropriate fields to floats
    studyspecs[['NPRAM', 'ZSTAT', 'P']] = studyspecs[['NPRAM', 'ZSTAT', 'P']].astype(float)

    # turn dataframes into dictionaries to process in db insertion
    filters = geneids.to_dict(orient='records')
    updates = studyspecs.to_dict(orient='records')

    bulk = list()
    for query_filter, query_update in zip(filters, updates):
        # modify the update query so that it contains the study ID and the push operator
        query_update = {'META': query_update.update({'META': query_update})}
        # add the update query to the bulk pipeline
        bulk.append(UpdateOne(query_filter, query_update, upsert=True))

    bulk_exec(col, bulk)