__author__ = 'olga'


def three_prime_motif(coordinate, strand, into_exon, into_intron,
                      debug=False):
    """
    Returns 20 bases of end of the intron, and 3 bases of the beginning of the exon
    """
    # If both are None, return None
    if into_exon is None and into_intron is None:
        return None, None

    # If one is none, replace with 1
    into_exon = 0 if into_exon is None else into_exon
    into_intron = 0 if into_intron is None else into_intron
    if debug:
        print 'into_intron', into_intron
        print 'into_exon', into_exon
    if strand == '+':
        return coordinate - into_intron - 1, coordinate + into_exon - 1
    elif strand == '-':
        return coordinate - into_exon - 1, coordinate + into_intron - 1
    else:
        raise TypeError('%s is not a valid strand identifier' % strand)


def five_prime_motif(coordinate, strand, into_exon, into_intron, debug=False):
    """
    Returns 3 bases of exon, and 6 bases of beginning of intron
    """
    if into_exon is None and into_intron is None:
        return None, None
    into_exon = 0 if into_exon is None else into_exon
    into_intron = 0 if into_intron is None else into_intron
    if debug:
        print 'into_intron', into_intron
        print 'into_exon', into_exon
    if strand == '+':
        return coordinate - into_exon, coordinate + into_intron
    elif strand == '-':
        return coordinate - into_intron - 1, coordinate + into_exon - 1
    else:
        raise TypeError('%s is not a valid strand identifier' % strand)


def motif_bed_line(chrom, start, stop,
                   trio, strand, ordinal,
                   name, debug=False):
    """
    @param ordinal: string, 'first', 'second', or 'third'. which exon in the trio
    @param name: descriptive string, e.g. '5prime' or '3prime'
    """
    if debug:
        print 'start:', start, ' stop:', stop
    if start is not None and stop is not None:
        return '%s\t%d\t%d\t%s\t.\t%s\n' % (
            chrom, start, stop,
            '%s|%s|%s' % (trio, ordinal, name), strand)


def intron_motif(intron_start, intron_stop, strand):
    start = intron_start if strand == '+' else intron_stop
    stop = intron_stop if strand == '+' else intron_start
    return start, stop


def add_line(lines, line):
    if line is not None:
        lines.append(line)
    return lines


def extract_motifs(trios_to_exons, name,
                   first_donor_into_exon=None,
                   first_donor_into_intron=None,
                   second_acceptor_into_exon=None,
                   second_acceptor_into_intron=None,
                   second_donor_into_exon=None,
                   second_donor_into_intron=None,
                   third_acceptor_into_exon=None,
                   third_acceptor_into_intron=None,

                   # these last three are only true/false
                   alt_exon=False,
                   constitutive_intron=False,
                   upstream_alt_intron=False,
                   downstream_alt_intron=False,

                   debug=False):
    """
    given a dictionary of MISO exon trios, extract the locations of things that
    are some integer away from exons, into the introns. Returns a single
    bed-formatted strings with proper tabs and newlines and everything.

    FYI: for the MISO annotations, the exons in the trios are always in the
    order exon1,exon2,exon3 regardless of the strand.

    @param trios_to_exons: dict of a string of trios, to the individual exons
    split on the @ sign
    @param first_donor_into_exon:
    @param first_donor_into_intron:
    @param second_acceptor_into_exon:
    @param second_acceptor_into_intron:
    @param second_donor_into_exon:
    @param second_donor_into_intron:
    @param third_acceptor_into_exon:
    @param third_acceptor_into_intron:
    @param alt_exon:
    @param constitutive_intron:
    @param upstream_alt_intron:
    @param downstream_alt_intron:
    @param name is a descriptor which will be added to the bed line. e.g. if
    'name' is 'upstream1000', the full ID will be the trio + the name.

    For example:
    >>> trios_to_exons = dict([('chr3:53274267:53274364:-@chr3:53271813:53271836:-@chr3:53268999:53269190:-',
    ...   ['chr3:53274267:53274364:-',
    ...    'chr3:53271813:53271836:-',
    ...    'chr3:53268999:53269190:-']),
    ...  ('chr2:9002720:9002852:-@chr2:9002401:9002452:-@chr2:9000743:9000894:-',
    ...   ['chr2:9002720:9002852:-',
    ...    'chr2:9002401:9002452:-',
    ...    'chr2:9000743:9000894:-']),
    ...  ('chr1:160192441:160192571:-@chr1:160190249:160190481:-@chr1:160188639:160188758:-',
    ...   ['chr1:160192441:160192571:-',
    ...    'chr1:160190249:160190481:-',
    ...    'chr1:160188639:160188758:-']),
    ...  ('chr7:100473194:100473333:+@chr7:100478317:100478390:+@chr7:100478906:100479034:+',
    ...   ['chr7:100473194:100473333:+',
    ...    'chr7:100478317:100478390:+',
    ...    'chr7:100478906:100479034:+']),
    ...  ('chr4:55124924:55124984:+@chr4:55127262:55127579:+@chr4:55129834:55130094:+',
    ...   ['chr4:55124924:55124984:+',
    ...    'chr4:55127262:55127579:+',
    ...    'chr4:55129834:55130094:+'])])
    >>> upstream1000 = extract_motifs(trios_to_exons,name='intron|upstream1000',
    ...                                      second_acceptor_into_intron=1000)
    >>> for l in upstream1000.splitlines()[:4]:
    ...    print l #doctest: +NORMALIZE_WHITESPACE
    chr3	53271835	53272835	chr3:53274267:53274364:-@chr3:53271813:53271836:-@chr3:53268999:53269190:-|second|acceptor|intron|upstream1000	.	-
    chr2	9002451	9003451	chr2:9002720:9002852:-@chr2:9002401:9002452:-@chr2:9000743:9000894:-|second|acceptor|intron|upstream1000	.	-
    chr1	160190480	160191480	chr1:160192441:160192571:-@chr1:160190249:160190481:-@chr1:160188639:160188758:-|second|acceptor|intron|upstream1000	.	-
    chr7	100477316	100478316	chr7:100473194:100473333:+@chr7:100478317:100478390:+@chr7:100478906:100479034:+|second|acceptor|intron|upstream1000	.	+

    For example, if first_donor_5p=10, then we'll extract 10 bases into the
    intron, from the first exon. This takes care of positive/negative stuff
    correctly too.

          first_donor   second_acceptor     second_donor     third_acceptor
    [ exon1 ]------------>-----------[ exon2 ]---------->--------------[ exon3 ]

    or, on the negative strand:

          third_acceptor   second_donor     second_acceptor     first_donor
    [ exon3 ]------------<-----------[ exon2 ]----------<--------------[ exon1 ]
    """
    lines = []
    for trio, exons in sorted(trios_to_exons.iteritems()):
        motif_start, motif_stop = None, None
        first_intron_stop, second_intron_stop = None, None
        second_intron_start, third_intron_start = None, None

        for i, exon in enumerate(exons):
        #             print i
            pieces = exon.split(':')
            chrom = pieces[0]
            exon_start = int(pieces[1])
            exon_stop = int(pieces[2])
            strand = pieces[3]


            if i == 0:
                first_intron_stop = exon_stop if strand == '+' else exon_start

                motif_start, motif_stop = \
                    five_prime_motif(first_intron_stop, strand,
                                     first_donor_into_exon,
                                     first_donor_into_intron, debug=debug)

                lines = add_line(lines, motif_bed_line(chrom, motif_start,
                                                       motif_stop,
                                                       trio, strand,
                                                       'first|donor',
                                                       name))
            if i == 1:
                if alt_exon:
                    if exon_start != exon_stop:
                        #continue
                        #print 'exon_start', exon_start
                        #print 'exon_stop', exon_stop
                        stop = exon_stop if strand == '+' else exon_stop
                        start = exon_start - 1 if strand == '+' else \
                            exon_start - 1
                        lines.append('%s\t%s\t%s\t%s\t%s\t%s\n'
                                     % (
                            chrom, start, stop, trio, '.', strand))
                        continue

                second_intron_stop = exon_start if strand == '+' else exon_stop
                motif_start, motif_stop = \
                    three_prime_motif(second_intron_stop, strand,
                                      second_acceptor_into_exon,
                                      second_acceptor_into_intron, debug=debug)
                lines = add_line(lines, motif_bed_line(chrom, motif_start,
                                                       motif_stop,
                                                       trio, strand,
                                                       'second|acceptor',
                                                       name))

                second_intron_start = exon_stop if strand == '+' else exon_start
                motif_start, motif_stop = \
                    five_prime_motif(second_intron_start, strand,
                                     second_donor_into_exon,
                                     second_donor_into_intron, debug=debug)
                lines = add_line(lines, motif_bed_line(chrom, motif_start,
                                                       motif_stop,
                                                       trio, strand,
                                                       'second|donor',
                                                       name))

            if i == 2:
                third_intron_start = exon_start if strand == '+' else exon_stop
                motif_start, motif_stop = \
                    three_prime_motif(third_intron_start, strand,
                                      third_acceptor_into_exon,
                                      third_acceptor_into_intron, debug=debug)
                lines = add_line(lines, motif_bed_line(chrom, motif_start,
                                                       motif_stop,
                                                       trio, strand,
                                                       'third|acceptor',
                                                       name, debug=debug))
        if constitutive_intron:
            exon_start, exon_stop = intron_motif(first_intron_stop,
                                                 third_intron_start, strand)
            lines = add_line(lines, motif_bed_line(chrom, exon_start,
                                                   exon_stop + 1,
                                                   trio, strand,
                                                   'intron|constitutive',
                                                   name))
        if upstream_alt_intron:
            exon_start, exon_stop = intron_motif(first_intron_stop,
                                                 second_intron_start,
                                                 strand)
            lines = add_line(lines, motif_bed_line(chrom, exon_start,
                                                   exon_stop + 1,
                                                   trio, strand,
                                                   'intron|upstream',
                                                   name))
        if downstream_alt_intron:
            exon_start, exon_stop = intron_motif(second_intron_stop,
                                                 third_intron_start, strand)
            lines = add_line(lines, motif_bed_line(chrom, exon_start,
                                                   exon_stop + 1,
                                                   trio, strand,
                                                   'intron|downstream',
                                                   name))
    return ''.join(lines)