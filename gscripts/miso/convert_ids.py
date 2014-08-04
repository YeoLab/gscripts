from collections import defaultdict
import sys

import numpy as np
import gffutils
import pandas as pd


def miso_exon_to_gencode_exon(exon):
    """Convert a single miso exon to gencode gffutils database exon id

    >>> # Skipped exon (SE) or Mutually exclusive exon (MXE) ID
    >>> miso_exon_to_gencode_exon('chr2:9624561:9624679:+')
    'exon:chr2:9624561-9624679:+'
    >>> # Alt 5' splice site (pick first of alternative)
    >>> miso_exon_to_gencode_exon('chr15:42565276:42565087|42565161:-')
    'exon:chr15:42565276-42565087:-'
    >>> # Alt 3' splice site (pick first of alternative)
    >>> miso_exon_to_gencode_exon('chr2:130914199|130914248:130914158:-')
    'exon:chr2:130914199-130914158:-'
    >>> # Retained intron: start, stop separated by '-' instead of ':'
    >>> miso_exon_to_gencode_exon('chr1:906259-906386:+')
    'exon:chr1:906259-906386:+'
    """
    return 'exon:{}:{}-{}:{}'.format(*miso_exon_to_coords(exon))


def miso_id_to_exon_ids(miso_id):
    """Convert a MISO-style alternative event ID

    Split on the pipe ("|") to account for Alt 5'/3' splice site events

    # A skipped exon (SE) ID
    >>> miso_id_to_exon_ids('chr2:9624561:9624679:+@chr2:9627585:9627676:+@chr2:9628276:9628591:+')
    ['exon:chr2:9624561-9624679:+', 'exon:chr2:9627585-9627676:+', 'exon:chr2:9628276-9628591:+']
    >>> # A mutually exclusive (MXE) ID
    >>> miso_id_to_exon_ids('chr16:89288500:89288591:+@chr16:89289565:89289691:+@chr16:89291127:89291210:+@chr16:89291963:89292039:+')
    ['exon:chr16:89288500-89288591:+', 'exon:chr16:89289565-89289691:+', 'exon:chr16:89291127-89291210:+', 'exon:chr16:89291963-89292039:+']
    >>> # An Alt 5' splice site (A5SS) ID
    >>> miso_id_to_exon_ids("chr15:42565276:42565087|42565161:-@chr15:42564261:42564321:-")
    ['exon:chr15:42565276-42565087:-', 'exon:chr15:42564261-42564321:-']
    >>> # An Alt 3' splice site (A3SS) ID
    >>> miso_id_to_exon_ids('chr2:130914824:130914969:-@chr2:130914199|130914248:130914158:-')
    ['exon:chr2:130914824-130914969:-', 'exon:chr2:130914199-130914158:-']
    >>> # A retained intron (RI) ID
    >>> miso_id_to_exon_ids('chr1:906066-906138:+@chr1:906259-906386:+')
    ['exon:chr1:906066-906138:+', 'exon:chr1:906259-906386:+']
    """
    return map(miso_exon_to_gencode_exon, miso_id.split('@'))


def miso_exon_to_coords(exon):
    """Convert a miso exon to gffutils coordinates

    >>> miso_exon_to_coords('chr2:130914824:130914969:-')
    ('chr2', '130914824', '130914969', '-')
    >>> # Alt 5' SS - pick the first of the alternative ends
    >>> miso_exon_to_coords('chr15:42565276:42565087|42565161:-')
    ('chr15', '42565276', '42565087', '-')
    >>> # Alt 3' SS - pick the first of the alternative starts
    >>> miso_exon_to_coords('chr2:130914199|130914248:130914158:-')
    ('chr2', '130914199', '130914158', '-')
    >>> # Retained intron
    >>> miso_exon_to_coords('chr1:906066-906138:+')
    ('chr1', '906066', '906138', '+')
    """
    strand = exon[-1]
    coords = map(lambda x: x.split('|')[0],
                 exon.split(':'))
    if '-' in coords[1]:
        start, stop = coords[1].split('-')
        coords = coords[0], start, stop, strand
    return coords[0], coords[1], coords[2], strand


def convert_miso_ids_to_everything(miso_ids, db,
                                   event_type,
                                   out_dir):
    """Given a list of miso IDs and a gffutils database, pull out the
    ensembl/gencode/gene name/gene type/transcript names, and write files
    into the out directory. Does not return a value.

    Parameters
    ----------
    miso_ids : list of str
        Miso ids to convert
    db : gffutils.FeatureDB
        gffutils feature database created from a gtf file
    event_type : str
        The type of splicing event. This is used for naming only
    out_dir : str
        Where to write the files to.
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

    n_miso_ids = len(miso_ids)
    sys.stdout.write('Converting {} {} miso ids using {} gffutils database '
                     'into {}.'.format(n_miso_ids, event_type, str(db),
                                       out_dir))

    for i, miso_id in enumerate(miso_ids):
        if i % 100 == 0:
            sys.stdout.write('On {}/{} {} miso ids'.format(i, n_miso_ids,
                                                           event_type))

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
                    map(lambda x: x.split('.')[0],
                        exon.attributes['gene_id']))
                gene_name.update(exon.attributes['gene_name'])
                gene_type.update(exon.attributes['gene_type'])
                gencode_transcript.update(exon.attributes['transcript_id'])
                ensembl_transcript.update(
                    map(lambda x: x.split('.')[0], exon.attributes[
                        'transcript_id']))
            except gffutils.FeatureNotFoundError:
                try:
                    #  not an exon, look for any overlapping transcripts here
                    prefix, chrom, startstop, strand = e.split(':')
                    start, stop = startstop.split('-')
                    transcripts = list(db.features_of_type('transcript',
                                                           strand=strand,
                                                           limit=(
                                                           chrom, int(start),
                                                           int(stop))))
                    # if there are overlapping transcripts...
                    if len(transcripts) == 0:
                        continue
                    else:
                        for transcript in transcripts:
                            gencode.update(transcript.attributes['gene_id'])
                            ensembl.update(
                                map(lambda x: x.split('.')[0],
                                    transcript.attributes['gene_id']))
                            gene_name.update(transcript.attributes['gene_name'])
                            gene_type.update(transcript.attributes['gene_type'])
                            gencode_transcript.update(
                                transcript.attributes['transcript_id'])
                            ensembl_transcript.update(
                                map(lambda x: x.split('.')[0],
                                    transcript.attributes['transcript_id']))
                except:
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
        sys.stdout.write('Wrote {}\n'.format(tsv))

    for name, d in to_misos.iteritems():
        tsv = '{}/{}_to_miso_{}.tsv'.format(out_dir, name, event_type)
        with open(tsv, 'w') as f:
            for k, v in d.iteritems():
                f.write('{}\t{}\n'.format(k, '\t'.join(v)))
        sys.stdout.write('Wrote {}\n'.format(tsv))
