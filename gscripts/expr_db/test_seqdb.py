import seq_db
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

engine = create_engine('mysql+pymysql://ppliu:some_pass@sauron.ucsd.edu/test')

seq_db.Base.metadata.create_all(engine)

Session = sessionmaker(bind=engine)
session = Session()
"""
try:
    s = seq_db.CLIPSeq()
    s.sample_name = 'test_sample'
    s.expr_name = 'test_experiment'
    s.file_location = 'some/where/on/the/server'
    s.species = 'hg19'
    s.collab = 'dr. sequencing'
    s.collab_institute = 'ucsd'

    session.add(s)
    session.commit()

except Exception as e: 
    
    print e
    session.rollback()
    session.commit()

"""
try:

    for expr in session.query(seq_db.SeqExpr).all():
        print expr.sample_name,
        print expr.project_name,
        print expr.check_file()

except Exception as e:
    print e
