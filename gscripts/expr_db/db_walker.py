"""
db_walker.py
a: ppliu

yeo lab 
ucsd 2013
"""

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey, \
    Boolean
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import select, and_, or_, not_
import pymysql
import argparse


class TableWalker(object):
    """
    """


    def __init__(self, table_walker_ID):
        self.table_walker_ID = table_walker_ID
        self.drop_lock = True
        self.engine = create_engine(
            'mysql+pymysql://ppliu:some_pass@sauron.ucsd.edu/{}' \
                .format(self.table_walker_ID))
        self.table_objects = {}
        self.metadata = MetaData()

        # Initialize Tables for SQLAlchemy here
        rnaseq = Table('rnaseq', self.metadata,
                       Column('exp_id', String(50), primary_key=True),
                       Column('species', String(50)),
                       Column('adapter', String(50)),
                       Column('compressed', Boolean()),
                       Column('paired', Boolean())
        )

        self.table_objects['rnaseq'] = rnaseq

        file_locations = Table('file_locations', self.metadata,
                               Column('exp_id', String(50), primary_key=True),
                               Column('file_path', String(100),
                                      primary_key=True),
                               Column('file_path2', String(100)),
        )

        self.table_objects['file_locations'] = file_locations
        return

    def addExpr(self, exp_id=None, species=None, adapter=None, compressed=False,
                paired=False, \
                file_path=None, file_path2=None):
        assert exp_id, "Specify experiment ID, keyword exp_id"
        assert species, "Specify species, keyword species"
        assert file_path, "Specify path to input, keyword file_path"

        file_ins = self.table_objects['file_locations'].insert().values(
            file_path=file_path, \
            file_path2=file_path2, exp_id=exp_id)

        param_ins = self.table_objects['rnaseq'].insert().values(exp_id=exp_id,
                                                                 species=species, \
                                                                 adapter=adapter,
                                                                 compressed=compressed,
                                                                 paired=paired)

        conn = self.engine.connect()
        conn.execute(file_ins)
        conn.execute(param_ins)
        return

    def unlock(self):
        """
        Function to remove the drop table lock on this class' table

        Args:
            None

        Returns:
            None
        """
        self.drop_lock = False

    def lock(self):
        """
        Function to apply the drop table lock on this class' table

        Args:
            None

        Returns:
            None
        """
        self.drop_lock = True

    def dropTables(self):
        """
        Drop this class' table from the database. Drop lock must be remove
        first with unlock() function.

        Args:
            None

        Returns:
            None
        """
        assert not self.drop_lock, "Cannot drop tables. Use unlock() first."
        if not self.drop_lock:
            for key in self.table_objects:
                self.table_objects[key].drop(self.engine)

        return

    def gatherParams(self, table, exp_id):

        conn = self.engine.connect()

        parameters = self.table_objects[table]
        param_dict = {}

        # select from parameters table
        s = select([parameters.c.parameter, \
                    parameters.c.value]). \
            where(parameters.c.exp_id == exp_id)

        for parameter, value in list(conn.execute(s)):
            param_dict[parameter] = value

        param_dict['exp_id'] = exp_id
        return param_dict

    def gatherFileLocs(self, exp_id):

        conn = self.engine.connect()
        # select from file locations table

        files_dict = {}

        files = self.table_objects['file_locations']

        s = select([files.c.file_path, files.c.file_path2]). \
            where(files.c.exp_id == exp_id)

        for file1, file2 in list(conn.execute(s)):
            files_dict['primary'] = file1
            files_dict['mate'] = file2

        files_dict['exp_id'] = exp_id
        return files_dict

    def createTables(self):

        try:

            self.metadata.create_all(self.engine)

        except Exception as e:
            print e
        return


