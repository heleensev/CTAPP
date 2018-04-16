#!/bin/bash
#$ -S /hpc/local/CentOS7/dhl_ec/software/sevpy/1/envs/py27/bin/python
#$ -cwd
#$ -l h_rt=03:45:00
#$ -l h_vmem=20G
#$ -o ~/logs/deploy/
#$ -e ~/logs/deploy/
#$ -N drugbank

import logging
import io
import re
import json
from bs4 import BeautifulSoup
from pymongo import ASCENDING


drugbank_file = "/hpc/local/CentOS7/dhl_ec/software/ctapp/data/drugbank_db.xml"
# drugbank_file = "/home/sevvy/Scripts/out/full_db/full_db_0.xml"

monconf_path = '/hpc/local/CentOS7/dhl_ec/software/sevpy/1/config/mon.conf'
logging.basicConfig(filename='/hpc/local/CentOS7/dhl_ec/software/ctapp/logs/drugbank.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger(__name__)


def get_mongo_client():
    from pymongo import MongoClient

    def get_host_info():
        with open(monconf_path, 'r') as mon_conf:
            lines = mon_conf.readlines()
        return lines[1].split('=')[1].strip()

    monip = get_host_info()
    client = MongoClient(host=monip)
    return client


def init_db(client, dbname, colname):
    db = client[dbname]
    collection = db[colname]
    return collection


def init_drugbank():
    pipeline = drugbank_parse()
    # launch mongo job
    # launch_mongo()
    # get client, create indexes
    # cli = get_mongo_client()
    # collection = init_db(cli, 'targetdb', 'drugbank')
    # collection.create_index([('targets.uniprot_id', ASCENDING)])
    #
    # bulk_insert(collection, pipeline, 'drugbank')


# write pipeline to database
def bulk_insert(collection, pipeline, name):
    try:
        logger.info('len pipe {}= {}'.format(name, len(pipeline)))
        if len(pipeline) > 1:
            collection.insert(pipeline)
    except Exception as e:
        logger.error(e)
        with open('drugbank.json', 'w') as out:
            json.dump(pipeline, out)


def drugbank_parse():
    bulk = list()
    drug = str()
    cnt = 0
    total = 0
    human = False
    org_pat = re.compile(r'(\s)*<organism>Human</organism>')

    drugbank = io.open(drugbank_file, 'r', encoding="UTF-8")

    # read the two header lines
    for i in range(2):
        line = drugbank.readline()

    # function for getting children tags from parent tags
    def get_attribute(parent, name, array=False):
        try:
            if array:
                return [x.text for x in parent.findAll(name)]
            else:
                return parent.find(name).text
        except AttributeError:
            pass

    line = True
    while line:
        line = drugbank.readline()
        if line == '</drug>\n':
            drug += line

            total += 1
            if human:
                drug_soup = BeautifulSoup(drug, 'lxml')
                target_array = list()

                for target in drug_soup.findAll('target'):
                    # print(target.findChildren(recursive=True))
                    if target.find('organism').text == 'Human':
                        target_doc = dict()
                        target_doc['name'] = get_attribute(target, 'name')
                        target_doc['actions'] = get_attribute(target, 'action', True)
                        target_doc['pubmed_ids'] = get_attribute(target, 'pubmed-id', True)
                        peptide = target.find('polypeptide')
                        target_doc['gene_name'] = get_attribute(peptide, 'gene-name')
                        target_doc['location'] = get_attribute(peptide, 'cellular-location')
                        target_doc['name'] = target.find('name').text
                        ids = target.findAll('external-identifier')
                        for id in ids:
                            if id.find('resource').text == 'UniProtKB':
                                target_doc['uniprot_id'] = id.find('identifier').text
                        target_array.append(target_doc)

                # if drug has human targets, gather drug info and add to the bulk pipeline
                if target_array:
                    drug_doc = dict()
                    drug_doc['drugbank_id'] = get_attribute(drug_soup, 'drugbank-id')
                    drug_doc['name'] = get_attribute(drug_soup, 'name')
                    drug_doc['description'] = get_attribute(drug_soup, 'description')
                    drug_doc['status'] = get_attribute(drug_soup, 'group')
                    drug_doc['indication'] = get_attribute(drug_soup, 'indication')
                    drug_doc['pharmacodynamics'] = get_attribute(drug_soup, 'pharmacodynamics')
                    drug_doc['mechanism'] = get_attribute(drug_soup, 'mechanism-of-action')
                    articles = drug_soup.find('articles')
                    drug_doc['pubmed_ids'] = get_attribute(articles, 'pubmed-id', True)
                    effects = drug_soup.find('snp-effects')
                    # snp_effects = get_attribute(effects, 'effect', True)
                    if effects:
                        snp_effects = effects.findAll('effect')
                        for snp in snp_effects:
                            try:
                                tag = snp.text.strip().split('\n')
                                symbol = tag[1]
                                rsid = tag[3]
                                allele = tag[4]
                                change = tag[5]
                                description = tag[6]
                                # make snp effect info a subdocument of corresponding target
                                for target in target_array:
                                    if target['gene_name'] == symbol:
                                        target['snp_effect'] = {'rs_id': rsid,
                                                                'allele': allele,
                                                                'change': change,
                                                                'description': description}
                            except ValueError as e:
                                logger.info(e)
                                logger.info(tag)
                                logger.info(drug_doc['drugbank_id'])

                    # add target info to the document
                    drug_doc['targets'] = target_array
                    # add the document to the pipeline
                    bulk.append(drug_doc)
                    cnt += 1
                human = False

            # empty the drug string for the next loop
            drug = str()

            if (total % 1000) == 0:
                print('total done: %s' % str(total))
        # if the line does not equal the end of drug tag, append to the drug string
        else:
            drug += line
            if org_pat.match(line):
                human = True

    print('bulk length: %s' % str(len(bulk)))
    with open('/hpc/local/CentOS7/dhl_ec/software/ctapp/data/drugbank/drug_bank_%s.json' % str(cnt), 'w') as out:
        json.dump(bulk, out)
    logger.info('%s processed drugs' % str(cnt))
    # insert the bulk pipeline to the database
    return bulk

if __name__ == "__main__":
    init_drugbank()
"""human drugs = 5789 """