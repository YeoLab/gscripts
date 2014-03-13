import numpy as np
import gffutils
from collections import defaultdict
import pandas as pd
import sys


def miso_exon_to_gencode_exon(exon):
    return 'exon:{}:{}-{}:{}'.format(*miso_exon_to_coords(exon))


def miso_id_to_exon_ids(miso_id):
    return map(miso_exon_to_gencode_exon, miso_id.split('@'))


def miso_exon_to_coords(exon):
    return exon.split(':')


def convert_miso_ids_to_everything(miso_ids, db,
                                   event_type,
                                   out_dir):
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
    out_dir = out_dir.rstrip('/')
    event_type = event_type.lower()

    miso_to_ensembl = {}
    miso_to_gencode = {}
    miso_to_gene_name = {}
    miso_to_gene_type = {}
    miso_to_ensembl_transcript = {}
    miso_to_gencode_transcript = {}

    ensembl_to_miso = defaultdict(list)
    gencode_to_miso = defaultdict(list)
    gene_name_to_miso = defaultdict(list)

    for miso_id in miso_ids:
        exons = miso_id_to_exon_ids(miso_id)

        gencode = set([])
        ensembl = set([])
        gene_name = set([])
        gene_type = set([])
        gencode_transcript = set([])
        ensembl_transcript = set([])
        for e in exons:
            try:
                exon = db[e]

                gencode.update(exon.attributes['gene_id'])
                ensembl.update(
                    map(lambda x: x.split('.')[0], exon.attributes['gene_id']))
                gene_name.update(exon.attributes['gene_name'])
                gene_type.update(exon.attributes['gene_type'])
                gencode_transcript(exon.attributes['transcript_id'])
                ensembl_transcript.update(
                    map(lambda x: x.split('.')[0], exon.attributes['gene_id'])
                )
            except gffutils.FeatureNotFoundError:
                continue

        if len(gencode) > 0:

            for ens in ensembl:
                ensembl_to_miso[ens].append(miso_id)
            for g in gene_name:
                gene_name_to_miso[g].append(miso_id)
            for g in gencode:
                gencode_to_miso[g].append(miso_id)

            ensembl = ','.join(ensembl)

            gencode = ','.join(gencode)
            gene_name = ','.join(gene_name)
            gene_type = ','.join(gene_type)

            gencode_transcript = ','.join(gencode_transcript)
            ensembl_transcript = ','.join(ensembl_transcript)

            miso_to_gencode[miso_id] = gencode
            miso_to_ensembl[miso_id] = ensembl
            miso_to_gene_name[miso_id] = gene_name
            miso_to_gene_type[miso_id] = gene_type

            miso_to_gencode_transcript[miso_id] = gencode_transcript
            miso_to_ensembl_transcript[miso_id] = ensembl_transcript

        else:
            miso_to_gencode[miso_id] = np.nan
            miso_to_ensembl[miso_id] = np.nan
            miso_to_gene_name[miso_id] = np.nan
            miso_to_gene_type[miso_id] = np.nan

    miso_tos = {'gencode_gene': miso_to_gencode,
                'ensembl_gene': miso_to_ensembl,
                'gene_name': miso_to_gene_name,
                'gene_type': miso_to_gene_type,
                'gencode_transcript': miso_to_gencode_transcript,
                'ensembl_transcript': miso_to_ensembl_transcript}

    to_misos = {'ensembl_gene': ensembl_to_miso,
                'gene_name': gene_name_to_miso,
                'gencode_gene': gencode_to_miso}

    for name, d in miso_tos.iteritems():
        df = pd.DataFrame.from_dict(d, orient='index')
        df.index.name = 'event_name'
        df.columns = [name]
        tsv = '{}/miso_{}_to_{}.tsv'.format(out_dir, event_type, name)
        df.to_csv(tsv, sep='\t')
        sys.stdout.write('Wrote {}'.format(tsv))

    for name, d in to_misos.iteritems():
        tsv = '{}/{}_to_miso_{}.tsv'.format(out_dir, name, event_type)
        with open(tsv, 'w') as f:
            for k, v in d.iteritems():
                f.write('{}\t{}\n'.format(k, '\t'.join(v)))
        sys.stdout.write('Wrote {}'.format(tsv))
