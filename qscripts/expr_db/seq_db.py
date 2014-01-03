from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, String, Boolean, Integer, ForeignKey

engine = create_engine('mysql+pymysql://ppliu:some_pass@sauron.ucsd.edu/test')
Base = declarative_base()

class SeqExpr(Base):
    __tablename__ = 'seqexpr'
    id = Column(Integer, primary_key = True)
    sample_name = Column(String(100))
    file_location = Column(String(200))
    pair_location = Column(String(200))
    species = Column(String(5))
    output_directory = ''
    compressed = False
    type = Column(String(50))
    collab = Column(String(50))
    

    __mapper_args__ = {
        'polymorphic_identity':'seqexpr',
        'polymorphic_on':type
    }
        

class RNASeq(SeqExpr):
    __tablename__ = 'rnaseq'    
    id = Column(Integer, ForeignKey('seqexpr.id'), primary_key=True)
    strand = Column(String(50))
    
    __mapper_args__ = {
        'polymorphic_identity':'rnaseq'
    }

