import seq_db
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
engine = create_engine('mysql+pymysql://ppliu:some_pass@sauron.ucsd.edu/test')


Session = sessionmaker(bind=engine)
session = Session()

try:
    s1 = seq_db.CLIPSeq()
    s1.sample_name = 'fto'
    s1.tech_rep = 1
    s1.project_name = 'fto'
    s1.file_location = '/projects/ps-yeolab/seqdata/20131003_Julia_FTO_CLIP/FTO_CLIP_1.fastq.gz'
    s1.species = 'hg19'

    s1.tube_label = 'fto1'
    s1.target = 'fto'
    s1.cell_source = 'npc'
    s1.experiment_date = '06/00/2013'
    s1.seq_kit = 'mabe227'
    s1.size_selection_min = 100
    s1.size_selection_max = 500
    s1.library_prep_generator = 'julia_nussbacher'
    s1.read_length = 50
    s1.sequencing_location = 'singapore'
    s1.sequencing_platform = 'gaii'
    s1.sequencing_submission_date = '07/12/2013'
    s1.clip_type = 'clip'    

    session.add(s1)

    s2 = seq_db.CLIPSeq()
    s2.sample_name = 'fto'
    s2.tech_rep = 2
    s2.project_name = 'fto'
    s2.file_location = '/projects/ps-yeolab/seqdata/20131003_Julia_FTO_CLIP/FTO_CLIP_2.fastq.gz'
    s2.species = 'hg19'

    s2.tube_label = 'fto2'
    s2.target = 'fto'
    s2.cell_source = 'npc'
    s2.experiment_date = '06/00/2013'
    s2.seq_kit = 'mabe227'
    s2.size_selection_min = 100
    s2.size_selection_max = 500
    s2.library_prep_generator = 'julia_nussbacher'
    s2.read_length = 50
    s2.sequencing_location = 'singapore'
    s2.sequencing_platform = 'gaii'
    s2.sequencing_submission_date = '07/12/2013'
    s2.clip_type = 'clip'
    
    session.add(s2)

    s3 = seq_db.CLIPSeq()
    s3.sample_name = 'fto'
    s3.tech_rep = 3
    s3.project_name = 'fto'
    s3.file_location = '/projects/ps-yeolab/seqdata/20131003_Julia_FTO_CLIP/FTO_CLIP_3.fastq.gz'
    s3.species = 'hg19'

    s3.tube_label = 'fto3'
    s3.target = 'fto'
    s3.cell_source = 'npc'
    s3.experiment_date = '06/00/2013'
    s3.seq_kit = 'mabe227'
    s3.size_selection_min = 100
    s3.size_selection_max = 500
    s3.library_prep_generator = 'julia_nussbacher'
    s3.read_length = 50
    s3.sequencing_location = 'singapore'
    s3.sequencing_platform = 'gaii'
    s3.sequencing_submission_date = '07/12/2013'
    s3.clip_type = 'clip'    

    session.add(s3)
    session.commit()

except Exception as e: 
    
    print e
    session.rollback()
    session.commit()


    
