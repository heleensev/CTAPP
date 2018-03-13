import pandas as pd
from bs4 import BeautifulSoup

pharm_ids = "/home/sevvy/Documents/Stage/Datasets/DrugBank/pharmacologically_active.csv"
xml_db = "/home/sevvy/Scripts/out/full_db/full_db_0.xml"

def id_parse():

    ids = pd.read_csv(pharm_ids, sep=',', header=0, names=['drugid', 'name', 'symbol', 'genbank protid', 'genbank genid',
                                                           'uniprotid', 'uniprottitle', 'pdb', 'genecard', 'genatlas',
                                                           'hgnc', 'species', 'drugids'])

    print(ids.head(10))
    print('org size '+str(ids.shape[0]))
    # ids.dropna(axis=1, inplace=True)
    namenans = ids['name'].dropna()
    print(namenans.shape[0])
    symnans = ids['symbol']
    print(symnans.shape[0])
    hgnans = ids['hgnc'].dropna()
    print(hgnans.shape[0])
    uninans = ids['uniprotid'].dropna()
    print(uninans.shape[0])
    genbankknans = ids['genbank genid'].dropna()
    print(genbankknans.shape[0])

    # print(ids.head(10))
    # print(ids['species'].head(20))

    ids.dropna(subset=['species'], inplace=True)
    nonhuman = ids[ids["species"].str.contains("Human")].shape[0]
    # print(nonhuman)
    hgnas = ids[ids['hgnc'].isnull()]
    # print(hgnas)
    hgnas = hgnas.dropna(subset=['species'])
    hgnanshuman = hgnas[hgnas['species'].str.contains("Human")]
    # print(hgnanshuman)
    # print(ids.head(50))

    #

def xml_db_parse():

    soup = BeautifulSoup(open(xml_db, 'r'), 'lxml')

    print(soup.title.string)
    drugs = soup.findAll('drug')

    for drug in drugs:

        drug_bank_id = drug.find('drugbank-id').text
        drug_name = drug.find('name').text
        description = drug.find('description').text
        status = drug.find('group').text

        articles = drug.find('articles')

        refs = drug.find('general-references')
        pm_ids = refs.findAll('pubmed-id')
        pm_ids = [ID.text for ID in pm_ids]

        indication = drug.find('indication').text
        pharmacodynamics = drug.find('pharmacodynamics').text
        mechanism = drug.find('mechanism-of-action').text

        atc_codes = list()
        for code in atc_codes:
            print(code)
        # print(drug.find('atc-codes'))
        gocodes= drug.find('go-classifiers')
        # for code in gocodes:
            # print(code.find('category'))
            # print(code.children())
        # print(drug.find_all())
        break
    # print(go_class)
xml_db_parse()

def count_drug_ids():

    xml = open('/home/sevvy/Documents/Stage/Datasets/DrugBank/full_db.xml', 'r')
    ids = list()

    line = True
    while line:
        line = xml.readline()
        if line.find('primary="true"') != -1:
            ids.append(line)

    print(len(ids))

def count_uniprot_ids():

    csv = '/home/sevvy/Documents/Stage/Datasets/DrugBank/drug_links.csv'
    import pandas as pd
    df = pd.read_csv(csv, sep=',', header=0)

    ids = df[['DrugBank ID', 'UniProt ID']]

    id_map = dict()
    for x in ids.itertuples():
        id_map.update({x[1]: x[2]})

    return id_map
    # select UniprotID, DrugbankID,

