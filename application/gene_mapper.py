#!/bin/bash
#$ -S /hpc/local/CentOS7/dhl_ec/software/sevpy/1/bin/python
#$ -cwd
#$ -o ~/logs/deploy/
#$ -e ~/logs/deploy/
#$ -N genemapper

import logging
import os
import subprocess
import sys
from glob import glob

import fire
import pandas as pd

sys.path.append('/home/dhl_ec/hseverin/deploy')
from application.GWASParser import GWASio
from application.Annotator import gene_db

chunker_path = '/hpc/local/CentOS7/dhl_ec/software/ctapp/application/map_chunker.py'
gene_mapper = '/hpc/local/CentOS7/dhl_ec/software/ctapp/application/gene_mapper.py'
MAGMA = '/hpc/local/CentOS7/dhl_ec/software/magma/magma'
gene_loc_data = '/hpc/local/CentOS7/dhl_ec/software/magma/NCBI37.3/NCBI37.3.gene.loc'
ref_dir = '/hpc/local/CentOS7/dhl_ec/software/magma/g1000_{0}/g1000_{0}'

magma_io = '/hpc/local/CentOS7/dhl_ec/software/ctapp/data/magma/'

logging.basicConfig(filename='../logs/genemapper.log', filemode='a', level=logging.DEBUG)
logger = logging.getLogger(__name__)


class MapGenes:

    def gene_analysis(self, id, N=None, eth='eur'):
        """
        gene annotation and analysis with MAGMA for each GWAS study
        :return:
        """
        # first concatenate the csv chunks of the particular study
        GWASio.concat(id)
        # if the study has a constant study size else the N column is used
        if N:
            n_arg = 'N={}'.format(N)
        else:
            n_arg = 'ncol=4'

        ref_data = ref_dir.format(eth.lower())

        def annotate(snp_loc, pref):
            if os.path.isfile(snp_loc):
                # call shell script to do annotation with MAGMA
                subprocess.call(
                    '{0} --annotate --snp-loc {2} --gene-loc {3} --out {4}'
                    .format(MAGMA, id, snp_loc, gene_loc_data, magma_io+'{}_{}'.format(pref, id)), shell=True)
                return True
            else:
                logger.info('SNP loc file {} does not exist'.format(snp_loc))

        def analyze(snp_loc):
            if os.path.isfile(snp_loc) and os.path.isfile(magma_io+'rs_{}.genes.annot'.format(id)):
                # call shell script to do gene analysis with MAGMA
                subprocess.call(
                    '{0} --bfile {1} --pval {2} use=SNP,P {3} --gene-annot {4}.genes.annot --genes-only --out {4}'
                    .format(MAGMA, ref_data, snp_loc, n_arg, magma_io+'rs_{}'.format(id), magma_io+'{}'.format(id)), shell=True)
                if os.path.isfile(magma_io+'{}.genes.out'.format(id)):
                    return True
                logger.error('gene analysis was not successful, check magma log')
            else:
                logger.error('SNP loc file {} or annotation file does not exist'.format(snp_loc))

        for prefix in ['rs', 'no_rs']:
            form_loc_data = '{0}{1}_{2}_complete.csv'.format(magma_io, prefix, id)

            # annotate and perform gene analyse on the data if the SNPs are mapped to 1000genome
            annot_pass = annotate(form_loc_data, prefix)
            annal_pass = False
            print('prefix: {}: annot_pass: {}'.format(prefix, annot_pass))
            if prefix == 'rs' and annot_pass:
                annal_pass = analyze(form_loc_data)
            if annot_pass:
                self.__split_gene_annot(id)
                gene_db.gene_update(id, prefix, annal_pass)

    @staticmethod
    def meta_analysis_chunker():
        """
        meta analysis of chromosome chunks of all the individual GWAS studies
        :return:
        """
        chr_no = [_ for _ in range(1, 24)]

        def get_params():
            if chr < 8 or chr == 23:
                return '20G', '00:20:00'
            else:
                return '10G', '00:10:00'

        # submit jobs per chromosome number, tweak runtime parameters according to chromosome size
        for chr in range(1, 24):
            vmem, rt = get_params()
            # call shell script to do gene analysis with MAGMA
            subprocess.call(
                'qsub -hold_jid MAGMA_* -l h_vmem={0} -l h_rt={1} -N MAGMA_meta{2} {3} MAGMA_meta_analysis {2}'
                .format(vmem, rt, chr, gene_mapper), shell=True)

    @staticmethod
    def MAGMA_meta_analysis(chr):

        filenames = sorted(glob(os.path.join(magma_io, '*_chr{}.csv'.format(chr))))
        cohorts = ' '.join(['{}.genes.out'.format(f.strip('.csv')) for f in filenames])
        subprocess.call(
            '{0} --meta genes={1} --out meta_chr{2}'.format(MAGMA, cohorts, chr), shell=True)
        if os.path.isfile(magma_io+'meta_chr{0}.genes.out'.format(chr)):
            gene_db.meta_update(chr)
        else:
            logger.error('MAGMA meta analysis was not successful, check logs')

    @staticmethod
    def __split_gene_annot(studyID):

        """
        split the output of MAGMA after gene analysis on chromosome number, so to
        perform meta analysis on multiple studies in chunked per chromosome
        :param studyID: ID per GWAS study
        :return:
        """

        path = os.path.join(magma_io, 'rs_{}.genes.out'.format(studyID))
        try:
            data = pd.read_csv(path, delim_whitespace=True, header=0)
        except FileNotFoundError:
            logger.error("file: {} not found in split_gene_annot, inspect magma annotation log".format(path))
            return

        chr_grouped = data.groupby(['CHR'])
        groups = [chr_grouped.get_group(x) for x in chr_grouped.groups]
        chrs = [key for key in chr_grouped.groups]

        for chr, group in zip(chrs, groups):
            path = os.path.join(magma_io, '{}_chr{}.csv'.format(studyID, chr))
            group.to_csv(path, header=True, sep='\t', index=False)


if __name__ == '__main__':
    fire.Fire(MapGenes)
