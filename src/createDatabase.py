#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Ronja Rösner

This module creates the SQL library from extracted data.

"""

import pandas as pd
import sqlite3, os
if __name__ == "__main__":
	from GeDaMa.src.downloadSpeciesData import ResolveData
else:
	from .downloadSpeciesData import ResolveData


column_names = {
			"taxonomy": [
				"IDX",
				"Kingdom",
				"Phylum",
				"Class",
				"taxOrder",
				"Family",
				"Genus",
				"Species",
				#"Subspecies",
				"ScientificName",
				"Authority",
				"Vernacular_Eng",
				"Vernacular_Ger"
			],
			"traits": [
				"IDX",
				"isMarine",
				"isBrackish",
				"isFresh",
				"isTerrestrial",
				"isAllWater",
				"isMarineFresh",
				"isExtinct"
			],
			"ids": [
				"IDX",
				"AccessionNumber",
				"usageKey",
				"IRMNG_ID",
				"AphiaID",
				"PESI_GUID"
			]
		}

def infer_sqlite_type(dtype):
	"""
	Infers the SQLite column type from a pandas dtype.

	Args:
		dtype (pandas dtype): The pandas dtype to infer from.

	Returns:
		str: The corresponding SQLite column type.
	"""
	if pd.api.types.is_integer_dtype(dtype):
		return "INTEGER"
	elif pd.api.types.is_numeric_dtype(dtype):
		return "NUMERIC"
	elif pd.api.types.is_float_dtype(dtype):
		return "REAL"
	elif pd.api.types.is_bool_dtype(dtype):
		return "BOOLEAN"
	elif pd.api.types.is_datetime64_any_dtype(dtype):
		return "DATETIME"
	elif pd.api.types.is_object_dtype(dtype):
		return "BLOB"
	else:
		return "TEXT"

# function used for removing nested
# lists in python using recursion
def remove_nestings(input_list: list, output_list: list):
	"""
	Gets a nested list as an input and returns a flattened version of that list.
	"""
	for i in input_list:
		if type(i) == list:
			remove_nestings(i, output_list)
		else:
			output_list.append(i)
	return(output_list)

class CreateDatabase():
	def __init__(self, value_dict: dict, db_file: str, log_function=print, stop_event=None):
		"""
		Class for creating a SQLite database from extracted data.

		Args:
			value_dict (dict): Dictionary containing a list of the following key-value pairs:
				- "sciNameList" (list): List of scientific names of species.
				- "accessionNumberList" (list): List of accession numbers for the species.
				- "getMissing" (bool): Flag to indicate whether to get missing accession numbers from NCBI database.
				- "append_data" (bool): Flag to indicate whether to append data to the existing database.
			log_function (function): Function for logging messages. Default is print.
		"""

		self.db_file = db_file
		self.log_function = log_function
		self.stop_event = stop_event
		self.append_db_file = value_dict["append_data"]
		self.getMissing = value_dict["getMissing"]
		self.sciNameList = value_dict["sciNameList"]
		self.sciNameCount = len(self.sciNameList)
		self.idx_list = range(0,self.sciNameCount)
		
		# get the list of accession numbers. if it is shorter than the list of species, return a list of only None items
		self.accNumberList = value_dict["accessionNumberList"] if len(value_dict["accessionNumberList"]) == self.sciNameCount else [None] * self.sciNameCount

	
	def makeDataframes(self):
		"""
		Function for making dataframes from the downloaded data.
		Creates three dataframes: taxonomy, traits, and ids.
		Each dataframe contains the data for each species in the list.

		Returns:
			dict: Dictionary containing the generated dataframes.
		"""
		def _get_data(sciName: str, accessionNumber: str):
			search_dict = {
				"sciName": sciName,
				"accessionNumber": accessionNumber,
				"getMissing": self.getMissing
			}
			data = ResolveData(search_dict, log_function=self.log_function)
			data_dict = data.combineDictionaries()
			return data_dict
		
		taxonomy_table = pd.DataFrame(columns=column_names["taxonomy"], index=self.idx_list)
		traits_table = pd.DataFrame(columns=column_names["traits"], index=self.idx_list)
		ids_table = pd.DataFrame(columns=column_names["ids"], index=self.idx_list)

		# if data should be appended to an existing database, check if the file exists and get the last IDX
		if os.path.isfile(self.db_file) and self.append_db_file:
			db_conn = sqlite3.connect(self.db_file)
			c = db_conn.cursor()
			c.execute("SELECT IDX FROM taxonomy")
			try:
				IDX_mod = max(c.fetchall())[0]+1
			except ValueError:
				IDX_mod = 1
			db_conn.close()
		else:
			IDX_mod = 0

		self.log_function(f"Downloading data for {self.sciNameCount} species...")
		for IDX, sciName in enumerate(self.sciNameList):
			# set the database index to the current index + the last index
			modified_IDX = IDX+IDX_mod
			if self.stop_event and self.stop_event.is_set():
				self.log_function("--- Download cancelled. ---\n")
				break
			
			self.log_function(f"[{IDX+1}]",end=' ')
			acc_number = self.accNumberList[IDX]
			data_dict = _get_data(sciName, acc_number)
			
			taxonomy_data = [
				modified_IDX,
				list(data_dict["Taxonomy"].values()),
				list(data_dict["Vernaculars"].values())
			]
			taxonomy_data = remove_nestings(taxonomy_data, [])
			taxonomy_table.iloc[IDX] = taxonomy_data

			traits_data = [
				modified_IDX,
				list(data_dict["Traits"].values())
			]
			traits_data = remove_nestings(traits_data, [])
			traits_table.iloc[IDX] = traits_data

			ids_data = [
				modified_IDX,
				list(data_dict["IDs"].values())
			]
			ids_data = remove_nestings(ids_data, [])
			ids_table.iloc[IDX] = ids_data

		out_dict = {
			"taxonomy": taxonomy_table,
			"traits": traits_table,
			"ids": ids_table
		}
		return out_dict
	
	def makeSQL(self):
		"""
		Function for creating a SQLite database from the dataframes.
		"""

		table_dict = self.makeDataframes()

		# remove database file if toggle is set to True, else print different statements
		if os.path.isfile(self.db_file):
			if self.append_db_file:
				self.log_function(f'\nDatabase already exists. Connecting...\n')
			else:
				os.remove(self.db_file)
				self.log_function(f'\nDatabase already exists. Removing and creating new file...\n')
		else:
			self.log_function(f'\nDatabase does not exist yet. Creating new file...\n')
		
		# establishes connection to the database and creating an empty file if there is none
		db_conn = sqlite3.connect(self.db_file)
		# creates new cursor object to interact with the database
		c = db_conn.cursor()

		for table_name, df in table_dict.items():
			columns = ', '.join(f'"{col}" {infer_sqlite_type(dtype)}' for col, dtype in df.dtypes.items())
			table_query = f"CREATE TABLE IF NOT EXISTS \"{table_name}\" ({columns})"
			c.execute(table_query)

			df.to_sql(table_name, db_conn, if_exists='append', index=False)

		db_conn.commit()
		db_conn.close()

		self.log_function(f"Data successfully loaded into internal database.")

def createNewDatabase(db_file: str, log_function=print):
	empty_value_dict = {
		"sciNameList": [],
		"accessionNumberList": [],
		"getMissing": False,
		"append_data": False
		}
	db_conn = sqlite3.connect(db_file)
	c = db_conn.cursor()
	c.execute("SELECT name FROM sqlite_master")
	table_names = c.fetchall()
	db_conn.close()

	if not table_names:
		print("Database is empty.")
		sql = CreateDatabase(empty_value_dict, db_file, log_function)
		sql.makeSQL()


if __name__ == "__main__":

	testValues = {
		"sciNameList": [
			"Oncorhynchus mykiss",
			"Calidris alpina",
			"Anguilla anguilla"
		],
		"accessionNumberList": [
			"GCF_002163495.1",
			"",
			"GCA_000695075.1"
		]
	}
	
	print(f"Test Species: {testValues['sciNameList']}")
	
	PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	TEST_DB_FILE = f"{PARENT_DIR}/data/test_database.db"

	createNewDatabase(TEST_DB_FILE)
	