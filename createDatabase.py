#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Ronja Rösner

This module creates the SQL library from extracted data.

"""

import pandas as pd
import sqlite3, os
from downloadSpeciesData import ResolveData, internetConnection

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DB_FILE = f"{SCRIPT_DIR}/data/species_database.db"

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

def get_data(sciName: str, accessionNumber: str = ""):
	data = ResolveData(sciName, accessionNumber)
	data_dict = data.combineDictionaries()
	return data_dict

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
	def __init__(self, value_dict: dict, log_function = print):
		self.log_function = log_function
		self.sciNameList = value_dict["sciNameList"]
		self.sciNameCount = len(self.sciNameList)
		self.idx_list = range(0,self.sciNameCount)
		
		# get the list of accession numbers. if it is shorter than the list of species, return a list of only None items
		self.accNumberList = value_dict["accessionNumberList"] if len(value_dict["accessionNumberList"]) == self.sciNameCount else [None] * self.sciNameCount

		
	
	def makeDataframes(self):
		taxonomy_table = pd.DataFrame(columns=column_names["taxonomy"], index=self.idx_list)
		traits_table = pd.DataFrame(columns=column_names["traits"], index=self.idx_list)
		ids_table = pd.DataFrame(columns=column_names["ids"], index=self.idx_list)

		self.log_function(f"Downloading data for {self.sciNameCount} species...")
		for IDX, sciName in enumerate(self.sciNameList):
			self.log_function(f"[{IDX+1}] Downloading data for {sciName}...")

			data_dict = get_data(sciName, self.accNumberList[IDX])
			
			taxonomy_data = [
				IDX,
				list(data_dict["Taxonomy"].values()),
				list(data_dict["Vernaculars"].values())
			]
			taxonomy_data = remove_nestings(taxonomy_data, [])
			taxonomy_table.iloc[IDX] = taxonomy_data

			traits_data = [
				IDX,
				list(data_dict["Traits"].values())
			]
			traits_data = remove_nestings(traits_data, [])
			traits_table.iloc[IDX] = traits_data

			ids_data = [
				IDX,
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
	
	def makeSQL(self, append_db_file: bool = True):
		table_dict = self.makeDataframes()

		# remove database file if toggle is set to True, else print different statements
		if os.path.isfile(DB_FILE):
			if append_db_file:
				self.log_function(f'\nDatabase already exists. Connecting...\n')
				
			else:
				os.remove(DB_FILE)
				self.log_function(f'\nDatabase already exists. Removing and creating new file...\n')
		else:
			self.log_function(f'\nDatabase does not exist yet. Creating new file...\n')
		
		# establishes connection to the database and creating an empty file if there is none
		db_conn = sqlite3.connect(DB_FILE)
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

	sql = CreateDatabase(testValues)
	sql.makeSQL()