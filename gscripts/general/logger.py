"""
logger.py
a: ppliu

yeo lab
ucsd 2013
"""

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey,\
	Boolean, Float
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import select, and_, or_, not_
import pymysql

class logger(object):
	"""
	Logger class will act as decorator, meant to decorate functions
	that return a dictionary of values to write to database. Parameters that
	uniquely identify the table, database and primary key of the set of values
	to be inserted must be provided.
	"""

	def __init__(self, unique_id, table, database, write=False):

		self.unique_id = unique_id
		self.table = table
		self.database = database
		self.write = write

	def __call__(self, func):

		if not self.write:
			return func

		engine = create_engine('mysql+pymysql://user:password@sauron.ucsd.edu/{}'\
				.format(self.database))
		
		metadata = MetaData()

		# reflection on passed database to gather existing tables
		metadata.reflect(engine)

		def write():

			# create new table if not present in database
			if not self.table in metadata:
				to_store = func()
					
				new_columns = []
				new_columns.append(Column('id', String(50), primary_key=True))

				for key in to_store:
					if isinstance(to_store[key], basestring):
						new_columns.append(Column(key, String(50)))

					elif isinstance(to_store[key], int):
						new_columns.append(Column(key, Integer))

					elif isinstance(to_store[key], float):
						new_columns.append(Column(key, Float))

				new_table = Table(self.table, metadata, *new_columns)
				new_table.create(engine, checkfirst=True)

				to_store['id'] = self.unique_id
				conn = engine.connect()
				print to_store
				conn.execute(new_table.insert(to_store))

			# write to exisiting table
			else:
				try:
					write_table = Table(self.table, metadata, \
						autoload=True, \
						autoload_with=engine)

					to_store = func()
					to_store['id'] = self.unique_id
					conn = engine.connect()
					print to_store
					conn.execute(write_table.insert(to_store))

				except Exception as e:
					print e
					return func

		return write

# test case
some_args = [10223, 'some_tool', 'test_walker', True]
@logger(*some_args)
def foo():
	print 'regular func'
	return {'x':1, 'y':2000, 'z':2}

foo()
