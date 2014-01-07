from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, String, Boolean, Integer, ForeignKey
from collections import defaultdict
import itertools
import string
import random

engine = create_engine('mysql+pymysql://ppliu:some_pass@sauron.ucsd.edu/test')
Base = declarative_base()

class SeqExpr(Base):

    __tablename__ = 'seqexpr'
    id = Column(Integer, primary_key = True)
    sample_name = Column(String(100))
    expr_name = Column(String(100))
    file_location = Column(String(200))
    pair_location = Column(String(200))
    species = Column(String(50))
    type = Column(String(50))
    collab = Column(String(50))
    collab_institute = Column(String(50))

    tube_label = Column(String(50))
    target = Column(String(50))
    cell_source = Column(String(50))
    experiment_date = Column(String(50))
    seq_kit = Column(String(50))
    size_selection_min = Column(Integer())
    size_selection_max = Column(Integer())
    library_prep_generator = Column(String(50))
    library_prep_date = Column(String(50))
    read_length = Column(Integer())
    sequence_location = Column(String(50))
    sequence_platform = Column(String(50))

    compressed = False

    __mapper_args__ = {
        'polymorphic_identity':'seqexpr',
        'polymorphic_on':type
    }
        

class RNASeq(SeqExpr):
    
    __tablename__ = 'rnaseq'  
    id = Column(Integer, ForeignKey('seqexpr.id'), primary_key=True)
    strand = Column(String(50))
    adapter = Column(String(50))   

    __mapper_args__ = {
        'polymorphic_identity':'rnaseq'
    }

    def get_params(self):
        return (self.__tablename__, self.species, self.strand, self.adapter)


class CLIPSeq(SeqExpr):
    
    __tablename__ = 'clipseq'
    id = Column(Integer, ForeignKey('seqexpr.id'), primary_key=True)

    __mapper_args__ = {
        'polymorphic_identity':'clipseq'
    }
    
    def get_params(self):
        return (self.__tablename__, self.species)


    
def manifest_setup(expr_list, working_dir='/home/ppliu/scratch/', filename='manifest.txt'):

    f = open(working_dir+filename, 'w')

    for expr in expr_list:
    
        if expr.__tablename__ == 'clipseq':
        
            f.write(expr.file_location)
            f.write('\t')
            f.write(expr.species)
            f.write('\n')

    f.close()

    pass
        


def q_setup(expr_list, working_dir='~/scratch/', session_name=None, queue='home', group='yeo'):

    star_genome = {
        'hg19': '/projects/ps-yeolab/genomes/hg19/star/',
        'mm9': '/projects/ps-yeolab/genomes/mm9/star/',
        'ce10': '/projects/ps-yeolab/genomes/ce10/star/',
        'dm3': '/projects/ps-yeolab/genomes/dm3/star/'
    }
    
    tags_annotation = {
        'hg19': '/projects/ps-yeolab/genomes/hg19/gencode_v17/gencode.v17.annotation.exons.bed',
        'mm9': '/projects/ps-yeolab/genomes/mm9/Mus_musculus.NCBIM37.64.fixed.exons.bed'
    }
    
    chr_sizes = {
        'hg19': '/projects/ps-yeolab/genomes/hg19/hg19.chrom.sizes',
        'mm9': '/projects/ps-yeolab/genomes/mm9/mm9.chrom.sizes'
    }

    scala_script = {
        'rnaseq': '~/gscripts/qscripts/analyze_rna_seq.scala',
        'clipseq': '~/gscripts/qscripts/analyze_clip_seq.scala'
    }


    params = [[expr.get_params(), expr] for expr in expr_list]
    for k, v in itertools.groupby(params, lambda x:"".join(x[0])):


        if session_name == None:

            sh_name = ''
            filelist_name = ''
            filelist_name = random.choice(string.letters)
            filelist_name += random.choice(string.letters)
            filelist_name += random.choice(string.letters)
            filelist_name += random.choice(string.letters)
            filelist_name += random.choice(string.letters)
            sh_name = filelist_name + '.sh'
            filelist_name += '.txt'
       
        else:
            sh_name = session_name + '.sh'
            filelist_name = session_name + '.txt'
 
        filelist_name = working_dir + '/' +filelist_name
        sh_name = working_dir + '/' + sh_name
        print filelist_name
       
        fl = open(filelist_name, 'w') 
        for item in list(v):

            if item[1].__tablename__ == 'rnaseq':
                fl.write(item[1].file_location)
                if item[1].pair_location:
                    fl.write('\t')
                    fl.write(item[1].pair_location)
                fl.write('\n')

            elif item[1].__tablename__ == 'clipseq':
                fl.write(item[1].file_location)
                
                    
            
        fl.close()
        
        if expr.__tablename__ == 'rnaseq':
            command = 'java -Xms512m -Xmx512m -jar ~/dash/gatk/dist/Queue.jar'
            command += ' -S '+ scala_script[expr.__tablename__]
            command += ' --input ' + filelist_name
            command += ' --species ' + expr.species
            command += ' --chr_sizes ' + chr_sizes[expr.species]
            command += ' --tags_annotation ' + tags_annotation[expr.species]
            command += ' --star_genome_location ' + star_genome[expr.species]
            command += ' --adapter ' + expr.adapter
            command += ' -jobQueue ' + queue
            command += ' -jobNative ' + '\"group_list={}-group\"'.format(group)
            command += ' --compressed '
            command += ' -run '
            command += ' -qsub '

        elif expr.__tablename__ == 'clipseq':
            
            command = 'java -Xms512m -Xmx512m -jar ~/dash/gatk/dist/Queue.jar'
            command += ' -S '+ scala_script[expr.__tablename__]
            command += ' --input ' + filelist_name
            command += ' --species ' + expr.species
            command += ' --chr_sizes ' + chr_sizes[expr.species]
            command += ' --tags_annotation ' + tags_annotation[expr.species]
            command += ' --star_genome_location ' + star_genome[expr.species]
            command += ' --adapter ' + 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
            command += ' --adapter ' + 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
            command += ' --adapter ' + 'ATCTCGTATGCCGTCTTCTGCTTG'
            command += ' --adapter ' + 'CGACAGGTTCAGAGTTCTACAGTCCGACGATC'
            command += ' --adapter ' + 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
            command += ' -jobQueue ' + queue
            command += ' -jobNative ' + '\"group_list={}-group\"'.format(group)
            command += ' --compressed '
            command += ' -run '
            command += ' -qsub '

        sh = open(sh_name, 'w')
        print k
        print command
        sh.write(command)
        sh.close()


   
