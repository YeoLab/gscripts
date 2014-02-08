import numpy as np
import gffutils
from collections import defaultdict, Counter


def miso_exon_to_gencode_exon(exon):
    return 'exon:{}:{}-{}:{}'.format(*miso_exon_to_coords(exon))


def miso_id_to_exon_ids(miso_id):
    return map(miso_exon_to_gencode_exon, miso_id.split('@'))


def miso_exon_to_coords(exon):
    return exon.split(':')


def convert_miso_id_to_everything(miso_ids, db):
    """
    Given a list of miso IDs and a gffutils database, pull out the
    ensembl/gencode/gene name/gene type/transcript names

    @param miso_ids:
    @type miso_ids:
    @param db:
    @type db:
    @return:
    @rtype:
    """

    miso_to_ensembl = {}
    miso_to_gencode = {}
    miso_to_gene_name = {}
    miso_to_gene_type = {}

    ensembl_to_miso = defaultdict(list)
    gene_name_to_miso = defaultdict(list)

    for miso_id in miso_ids:
        exons = miso_id_to_exon_ids(miso_id)

        gencode = set([])
        ensembl = set([])
        gene_name = set([])
        gene_type = set([])
        for e in exons:
            try:
                exon = db[e]

                gencode.update(exon.attributes['gene_id'])
                ensembl.update(
                    map(lambda x: x.split('.')[0], exon.attributes['gene_id']))
                gene_name.update(exon.attributes['gene_name'])
                gene_type.update(exon.attributes['gene_type'])
            except gffutils.FeatureNotFoundError:
                continue

        if len(gencode) > 0:

            ensembl = ','.join(ensembl)

            gencode = ','.join(gencode)
            gene_name = ','.join(gene_name)
            gene_type = ','.join(gene_type)

            miso_to_gencode[miso_id] = gencode
            miso_to_ensembl[miso_id] = ensembl
            miso_to_gene_name[miso_id] = gene_name
            miso_to_gene_type[miso_id] = gene_type

            ensembl_to_miso[ensembl].append(miso_id)
            gene_name_to_miso[gene_name].append(miso_id)
        else:
        #         print gencode.most_common(1)[0][0]
            print miso_id
            miso_to_gencode[miso_id] = np.nan
            miso_to_ensembl[miso_id] = np.nan
            miso_to_gene_name[miso_id] = np.nan
            miso_to_gene_type[miso_id] = np.nan

