__author__ = 'lovci'


"""a tool to extract ranges from indexed .maf file, adapted from bx-python / scripts / maf_tile.py and \
http://biopython.org/wiki/Phylo_cookbook"""

import bx
import bx.align.maf
from Bio import Phylo
import pyfasta
from collections import defaultdict
import os, sys
import socket

if "tscc" in socket.gethostname():
    conservation_basedir = "/projects/ps-yeolab/conservation"
    genome_basedir = "/projects/ps-yeolab/genomes"
else:
    #haven't set up this system yet... download and index .maf files
    raise NotImplementedError

class MafRangeGetter(object):

    def intervals_from_mask(self, mask ):
        start = 0
        last = mask[0]
        for i in range( 1, len( mask ) ):
            if mask[i] != last:
                yield start, i, last
                start = i
                last = mask[i]
        yield start, len(mask), last

    def tile_interval(self, chrom, start, end, strand):

        """where the magic happens... tile blocks from a .maf file to make a contiguous alignment
        for species other than the pivot species, reference (chromosome) names are masked, as the blocks
        from which a tiled interval may originate may not be on the same reference in an ortholog"""

        base_len = end - start
        ref_src = self.species + "." + chrom
        blocks = self.index.get( ref_src, start, end )
        # From low to high score
        blocks.sort( lambda a, b: cmp( a.score, b.score ) )

        mask = [ -1 ] * base_len

        ref_src_size = None
        if len(blocks) > 0 and blocks[0] is not None:
            for i, block in enumerate( blocks ):
                ref = block.get_component_by_src_start( ref_src )
                ref_src_size = ref.src_size
                assert ref.strand == "+"
                slice_start = max( start, ref.start )
                slice_end = min( end, ref.end )
                for j in range( slice_start, slice_end ):
                    mask[j-start] = i

        tiled = []
        for i in range( len( self.sources ) ): tiled.append( [] )
        for ss, ee, index in self.intervals_from_mask( mask ):
            if index < 0: #masked portions
                fetch_start = (start+ss-1)
                fetch_stop = (start +ee-1)
                x = self.fasta[chrom][fetch_start:fetch_stop]
                tiled[0].append(x)
                for row in tiled[1:]: #unaligned other species
                    row.append( "-" * ( ee - ss ) )
            else:
                slice_start = start + ss
                slice_end = start + ee
                block = blocks[index]
                ref = block.get_component_by_src_start( ref_src )
                sliced = block.slice_by_component( ref, slice_start, slice_end )
                sliced = sliced.limit_to_species( self.sources )
                sliced.remove_all_gap_columns()
                for i, src in enumerate( self.sources ):
                    comp = sliced.get_component_by_src_start( src )
                    if comp:
                        tiled[i].append( comp.text )
                    else:
                        tiled[i].append( "-" * sliced.text_size )
        aln = bx.align.Alignment()
        s = list()


        for i, name in enumerate( self.sources ):
            for n, s in enumerate(tiled[i]):
                tiled[i][n] = s.strip().replace("\s", "")
            text = str("".join( tiled[i]))
            size = len( text ) - text.count( "-" )
            if i == 0:
                if ref_src_size is None: ref_src_size = len(self.fasta[chrom])
                c = bx.align.Component( ref_src, start, end-start, "+", ref_src_size, text )
            else:
                c = bx.align.Component( name + ".fake", 0, size, "?", size, text )
            aln.add_component( c )
        if strand == "-":
            aln = aln.reverse_complement()
        return aln

    def _lookup_tree(self):
        """
        lookup dictionary for nodes
        from: http://biopython.org/wiki/Phylo_cookbook
        """
        names = {}
        for clade in self.tree.find_clades():
            if clade.name:
                if clade.name in names:
                    raise ValueError("Duplicate key: %s" % clade.name)
                names[clade.name] = clade
        self.treeNodes = names

    def _phylogenetic_distance_table(self):
        """precompute distances"""
        distance = defaultdict()
        for species in self.treeNodes.keys():
            bl = 0
            for b in self.tree.get_path(self.treeNodes[species]):

                bl += b.branch_length
            distance[species] = bl
        self.phyloD = distance
    def plot_tree(self):
        """show phylogenetic tree"""

        import pylab
        f = pylab.figure(figsize=(8,8))
        ax = f.add_subplot(111)
        y = Phylo.draw(self.tree, axes=ax)

class ce10MafRangeGetter(MafRangeGetter):

    def __init__(self):
        self.species = "ce10"
        super(MafRangeGetter,self).__init__()
        chrs = ["I", "II", "III", "IV", "V", "X", "M"]
        maf_files = [(os.path.join(conservation_basedir, "ce10_7way/", "multiz7way/", "chr" + c + ".maf.bz2"))\
                     for c in chrs]
        self.index = bx.align.maf.MultiIndexed( maf_files, keep_open=True, parse_e_rows=True, use_cache=True)
        self.fasta = pyfasta.Fasta(os.path.join(genome_basedir, "/ce10/chromosomes/all.fa"), flatten_inplace=True)
        treeFile = os.path.join(conservation_basedir, "ce10_7way/", "multiz7way/", "ce10.7way.nh")
        self.tree = Phylo.read(treeFile, 'newick')
        self.sources = [i.name for i in self.tree.get_terminals()]
        self.tree.root_with_outgroup({"name":"ce10"})
        self.tree.ladderize()
        self._lookup_tree() #make tree search-able
        self._phylogenetic_distance_table() #calculate the distance from each species to hg19


class hg19MafRangeGetter(MafRangeGetter):

    def __init__(self):
        self.species = "hg19"
        super(MafRangeGetter,self).__init__()
        chrs = map(str, range(1,23)) + ["X", "Y", "M"]
        maf_files = [os.path.join(conservation_basedir, "hg19_46way/zipped/chr" + c + ".maf.bz2")\
                     for c in chrs]
        self.index = bx.align.maf.MultiIndexed( maf_files, keep_open=True, parse_e_rows=True, use_cache=True )
        self.fasta = pyfasta.Fasta(os.path.join(genome_basedir, "hg19/chromosomes/all.fa"), flatten_inplace=True)

        treeFile = os.path.join(conservation_basedir, "hg19_46way/46way.corrected.nh")
        self.tree = Phylo.read(treeFile, 'newick')
        self.sources = [i.name for i in self.tree.get_terminals()]
        self.tree.root_with_outgroup({"name":"hg19"})
        self.tree.ladderize()

        self._lookup_tree() #make tree search-able
        self._phylogenetic_distance_table() #calculate the distance from each species to hg19

def revcom(s):
    import string
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV-', 'tgcayrkmvhdbTGCAYRKMVHDB-')
    s = s.translate(complements)[::-1]

    return s


def chop_maf(maf, maxSize=2500, overlap=6, id = "none", verbose=False):
    """make smaller bits (<=maxSize) from a big maf range with little bits overlapping by $overlap nt"""
    i = 0
    j = maxSize
    fullSize = len(maf.components[0].text)


    while i < fullSize:
        if i > 0 and verbose:
            pDone = 100*float(i)/fullSize
            sys.stderr.write( "chunking id:%s from %d-%d, %3.2f \r" %(id, i, j, pDone)  )
        yield maf.slice(i,j)

        i = j + -overlap
        j = i + maxSize


def test():
    """for gpratt"""
    getter = ce10MafRangeGetter()
    lin41 = getter.tile_interval("chrI", 9335957, 9341889, '-')
    assert lin41.slice(0,100).components[0].text == \
           'ATGGCGA---------CCATC------------GTGCCATGCTCATTGGAGAAAGAAGAAGGAGCAC---------CATCA-------------------'

