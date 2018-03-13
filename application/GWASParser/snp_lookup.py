# module for checking if SNPid is present in dbSNP/dbVar
import os
import requests
from numpy import NaN
from bs4 import BeautifulSoup as BS
from pyliftover import LiftOver

chainpath = '../data/chains/'
hg16hg19 = 'hg16ToHg19.over.chain.gz'
hg17hg19 = 'hg17ToHg19.over.chain.gz'
hg18hg19 = 'hg18ToHg19.over.chain.gz'
hg38hg19 = 'hg38ToHg19.over.chain.gz'

chains = {'hg16': hg16hg19, 'hg17': hg17hg19,
          'hg18': hg18hg19, 'hg38': hg38hg19}


def liftover_to_19(loc, build):
    floc = [loc.split(':')[0], loc.split(':')[1]]
    lo = LiftOver(os.path.join(chainpath, chains.get(build)))
    con_pos = lo.convert_coordinate(*floc)
    if con_pos:
        return int(con_pos[0][1])
    return NaN

def liftover_loc(study, loc):
    try:
        def liftover():
            lo = LiftOver(chain)
            con_pos = lo.convert_coordinate(*floc)[0]
            if con_pos:
                return con_pos[0][0], con_pos[0][1]

        floc = [loc.split(':')[0], loc.split(':')[1]]
        if study.build:
            chain = chains.get(study.build)
            if study.build != '38':
                floc = liftover()
            rs_id = query_loc(*floc)
            return rs_id
        else:
            return None
    except:
        # logger: not possible to liftover genomic position
        pass


def query_loc(chr, bp):

    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/' \
          'esearch.fcgi?db=snp&term={}[Base%20position]{}' \
          '[Chromosome]&term%22Homo+sapiens%22[Organism]'
    tries = 0
    while tries < 10:
        try:
            resp = requests.get(url.format(bp, chr))
            xml = resp.text
            soup = BS(xml, 'lxml')
            idlist = soup.find_all('idlist')
            hit = idlist[0].id
            if hit:
                return str(hit).split('<')[1].split('>')[1]
            break
        except requests.RequestException:
            tries += 1
            continue
        except:
            break


def query_snp(ids, db_name):
    import requests
    from bs4 import BeautifulSoup as bs
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={}&term={}[Accession]'

    for query in ids:
        go = True
        tries = 0
        while go and tries < 10:
            try:
                resp = requests.get(url.format(db_name, query))
                xml = resp.text
                soup = bs(xml, 'lxml')
                idlist = soup.find_all('idlist')
                hit = idlist[0].id
                if hit:
                    return query
            except requests.RequestException:
                tries += 1
            except:
                break