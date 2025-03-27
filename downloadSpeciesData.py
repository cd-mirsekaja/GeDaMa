#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Ronja Rösner

This module extracts data from multiple sources:
	Wikipedia
	the GBIF backbone
	the NCBI Genome database
"""

# import libraries
import os, requests, sqlite3, json
from pygbif import species as sp
import wikipediaapi as wiki
from Bio import Entrez
from time import sleep


def internetConnection():
	"""
	Function for checking if the user is connected to the internet.
	"""
	try:
		requests.get("https://api.gbif.org/", timeout=5)
		return True
	except requests.ConnectionError:
		return False

def restResponse(api_url: str):
		"""
		Function for saving a server response to a dictionary.
		Returns a dictionary of strings.
		"""

		response = requests.get(api_url)
		if response.status_code != 200:
			return None
		try:
			response_data = response.json()
			return response_data
		except json.JSONDecodeError:
			return None

def getAccessionNumbers(scientific_names, log_function=print, api_key=None, max_results=20):
	"""
	Retrieves the best genome accession (RefSeq preferred) from NCBI Assembly.
	This function was written with the help of Perplexity AI Sonar.
	
	Args:
		scientific_names (list): Organism scientific names
		api_key (str): NCBI API key for rate limits
		max_results (int): Max assemblies to check per species
	
	Returns:
		dict: {organism: [best_accession] or empty list}
	"""
	Entrez.email = "your_email@example.com"
	if api_key:
		Entrez.api_key = api_key
	
	accessions = {}
	log_function(f"Looking up Accessions for {len(scientific_names)} species...")
	for IDX, name in enumerate(scientific_names):
		log_function(f"[{IDX+1}] Searching for {name}...")
		try:
			# Get Taxonomy ID
			tax_handle = Entrez.esearch(db="taxonomy", term=name, retmax=1)
			tax_data = Entrez.read(tax_handle, validate=False)
			tax_handle.close()
			
			if not tax_data["IdList"]:
				accessions[name] = []
				continue
				
			tax_id = tax_data["IdList"][0]
			
			# Search Assembly database
			assembly_handle = Entrez.esearch(
				db="assembly",
				term=f"txid{tax_id}[Organism]",
				retmax=max_results,
				sort="date desc"  # Get newest first
			)
			assembly_data = Entrez.read(assembly_handle, validate=False)
			assembly_handle.close()
			
			assembly_ids = assembly_data["IdList"]
			if not assembly_ids:
				accessions[name] = []
				continue

			# Fetch assembly details
			best_acc = None
			for assembly_id in assembly_ids:
				summary_handle = Entrez.esummary(
					db="assembly", 
					id=assembly_id,
					retmode="xml"
				)
				summary = Entrez.read(summary_handle, validate=False)
				summary_handle.close()
				
				ds = summary["DocumentSummarySet"]["DocumentSummary"][0]
				primary_acc = str(ds["AssemblyAccession"])
				
				# Priority 1: RefSeq accession
				if primary_acc.startswith("GCF_"):
					best_acc = primary_acc
					break
				
				# Priority 2: Most recent GenBank
				if not best_acc and primary_acc.startswith("GCA_"):
					best_acc = primary_acc
				
				sleep(0.34)  # NCBI rate limit

			accessions[name] = best_acc if best_acc else ""
			if accessions[name] != "":
				log_function(f"Accession {best_acc} found!")
			else:
				log_function("No Accession Number found.")
			sleep(1)
			
		except Exception as e:
			print(f"Error processing {name}: {str(e)}")
			accessions[name] = []
	
	accessions_list = list(accessions.values())
	return accessions, accessions_list

def getScientificNames(acc_numbers: list, log_function=print):
	API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/"
	url_seg_accession = "genome/accession/"

	log_function(f"Looking up Scientific Names for {len(acc_numbers)} genomes...")
	sci_names = {}
	for IDX, acc_number in enumerate(acc_numbers):
		dataset_response = requests.get(API_URL+url_seg_accession+acc_number+"/dataset_report")
		if dataset_response.status_code == 200:
			sci_name = dataset_response.json()["reports"][0]["organism"]["organism_name"]
			sci_names[acc_number] = sci_name
			log_function(f"[{IDX+1}] {sci_name} found for {acc_number}")
		else:
			sci_names[acc_number] = ""
			log_function(f"[{IDX+1}] No name found for {acc_number}")
	
	sci_names_list = list(sci_names.values())
	return sci_names, sci_names_list
	

def getIdentification(sciName: str, accessionNumber: str):
	sciNameMod = '%20'.join(sciName.split())

	rest_url = "https://www.irmng.org/rest/IRMNG_IDByName/"
	api_url = (rest_url+sciNameMod)
	irmng_id = restResponse(api_url)

	rest_url = "https://www.eu-nomen.eu/portal/rest/PESIGUIDByName/"
	api_url = (rest_url+sciNameMod)
	guid_id = restResponse(api_url)

	rest_url = "https://www.marinespecies.org/rest/AphiaIDByName/"
	rest_mod = "?marine_only=false&extant_only=false"
	api_url = (rest_url+sciNameMod+rest_mod)
	aphia_id = restResponse(api_url)

	gbif_usageKey = sp.name_backbone(sciName)['usageKey']

	id_dict = {
		"AccessionNumber": accessionNumber,
		"GBIF-UsageKey": gbif_usageKey,
		"IRMNG-ID": irmng_id,
		"AphiaID": aphia_id,
		"PESI-GUID": guid_id
	}
	return id_dict

def getResponses(sciName: str, accessionNumber: str):
	"""
		Function for saving the server responses to all API calls. Calls to:
		- World Register of Marine Species (WORMS)
		- Interim Register of Marine and Nonmarine Genera (IRMNG)
		- Global Biodiversity Information Facility (GBIF)
		- Pan-European Species directories Infrastructure (PESI)
		- Wikipedia

		Returns a dictionary of dictionaries.
		"""
	
	sciNameMod = '%20'.join(sciName.split())
	id_dict = getIdentification(sciName, accessionNumber)

	# get GBIF response
	rest_url = "https://api.gbif.org/v1/species/"
	api_url = (f"{rest_url}{id_dict['GBIF-UsageKey']}")
	response = restResponse(api_url)
	try:
		gbif_response = response
	except TypeError:
		gbif_response = {}

	# get IRMNG response
	rest_url = "https://www.irmng.org/rest/AphiaRecordsByMatchNames?scientificnames%5B%5D="
	api_url = (rest_url+sciNameMod)
	response = restResponse(api_url)
	try:
		irmng_response = response[0][0]
	except TypeError:
		irmng_response = {}

	# get PESI response
	rest_url = "https://www.eu-nomen.eu/portal/rest/PESIRecordsByMatchTaxon/"
	api_url = (rest_url+sciNameMod)
	response = restResponse(api_url)
	try:
		pesi_response = response[0]
	except TypeError:
		pesi_response = {}

	# get WORMS response
	rest_url = "https://www.marinespecies.org/rest/AphiaRecordsByNames?scientificnames%5B%5D="
	rest_mod = "&like=false&marine_only=false&extant_only=false"
	api_url = (rest_url+sciNameMod+rest_mod)
	response = restResponse(api_url)
	try:
		worms_response = response[0][0]
	except TypeError:
		worms_response = {}
	
	# get GBIF vernaculars
	def _gbif_vernaculars():
		base_url="https://api.gbif.org/v1/species/"
		rest_service="/vernacularNames"
		api_url=(f"{base_url}{id_dict['GBIF-UsageKey']}{rest_service}")
		response = restResponse(api_url)
		try:
			gbif_ver_response = response
		except TypeError:
			gbif_ver_response = {}

		if gbif_ver_response!=None:
			engver_list = [result["vernacularName"] for result in gbif_ver_response["results"] if result["language"] == "eng"]
			gerver_list = [result["vernacularName"] for result in gbif_ver_response["results"] if result["language"] == "deu"]
		else:
			engver_list = []
			gerver_list = []

		if engver_list==[]:
			eng_out = ""
		else:
			eng_out = str(engver_list[0])
		
		if gerver_list == []:
			ger_out = ""
		else:
			ger_out = str(gerver_list[0])

		return {"English": eng_out, "German": ger_out}
	
	# get WORMS vernaculars
	def _worms_vernaculars():
		base_url = "https://www.marinespecies.org/rest/AphiaVernacularsByAphiaID/"
		api_url = (f"{base_url}{id_dict['AphiaID']}")
		response = restResponse(api_url)
		try:
			worms_ver_response = response
		except TypeError:
			worms_ver_response = {}

		if worms_ver_response!=None:
			engver_list = [entry["vernacular"] for entry in worms_ver_response if entry["language_code"] == "eng" or entry["language"] == "English"]
			gerver_list = [entry["vernacular"] for entry in worms_ver_response if entry["language_code"] == "deu" or entry["language"] == "German"]
		else:
			engver_list = []
			gerver_list = []
		
		if engver_list == []:
			eng_out = ""
		else:
			eng_out = str(engver_list[0])

		if gerver_list == []:
			ger_out = ""
		else:
			ger_out = str(gerver_list[0])

		return {"English": eng_out, "German": ger_out}
	
	def _wikipedia_vernaculars(sciName: str):
		wiki_en = wiki.Wikipedia('VernacularData','en')
		wiki_de = wiki.Wikipedia('VernacularData','de')

		if wiki_en.page(sciName).exists():
			pagetitle_en=wiki_en.page(sciName).displaytitle
			if "<i>" in pagetitle_en:
				pagetitle_en = ""
		else:
			pagetitle_en = ""

		if wiki_de.page(sciName).exists():
			pagetitle_de=wiki_de.page(sciName).displaytitle
			if "<i>" in pagetitle_de:
				pagetitle_de = ""
		else:
			pagetitle_de = ""
			
		return {"English": pagetitle_en, "German": pagetitle_de}

	# combine all responses into a dictionary and return it
	data_dict = {
		"WORMS": worms_response,
		"GBIF": gbif_response,
		"IRMNG": irmng_response,
		"PESI": pesi_response,
		"WORMS-Vernaculars": _worms_vernaculars(),
		"GBIF-Vernaculars": _gbif_vernaculars(),
		"Wikipedia-Vernaculars": _wikipedia_vernaculars(sciName)
	}
	return data_dict


class GetInformation:
	def __init__(self, sciName: str, accessionNumber: str):
		self.response_dict = getResponses(sciName, accessionNumber)
		self.id_dict = getIdentification(sciName, accessionNumber)
	
	# function for retrieving the taxonomic path from GBIF
	def gbifTaxpath(self):
		"""
		Function for getting the taxnonomic data from the GBIF response.
		Returns a dictionary of strings.
		"""

		gbif_data = self.response_dict["GBIF"]
		tkingdom, tphylum, tclass, torder, tfamily, tgenus, tspecies, tsciName, tauthority = "","","","","","","","",""
		if gbif_data != {}:
			if 'kingdom' in gbif_data:
				tkingdom = gbif_data['kingdom']
			if 'phylum' in gbif_data:
				tphylum = gbif_data['phylum']
			if 'class' in gbif_data:
				tclass = gbif_data['class']
			if 'order' in gbif_data:
				torder = gbif_data['order']
			if 'family' in gbif_data:
				tfamily = gbif_data['family']
			if 'genus' in gbif_data:
				tgenus = gbif_data['genus']
			if gbif_data["rank"] == "SPECIES":
				if 'scientificName' in gbif_data:
					tspecies = gbif_data['scientificName'].split(" ",2)[1]
					tauthority = gbif_data['scientificName'].split(" ",2)[2]
				if 'species' in gbif_data:
					tsciName = gbif_data['species']
			
		out_dict={
			"Kingdom": tkingdom,
			"Phylum": tphylum,
			"Class": tclass,
			"taxOrder": torder,
			"Family": tfamily,
			"Genus": tgenus,
			"Species": tspecies,
			"ScientificName": tsciName,
			"Authority": tauthority
			}
		return out_dict
	
	# function for retrieving the taxonomic path from IRMNG
	def irmngTaxpath(self):
		"""
		Function for getting the taxnonomic data from the IRMNG response.
		Returns a dictionary of strings.
		"""

		irmng_data = self.response_dict["IRMNG"]
		tkingdom, tphylum, tclass, torder, tfamily, tgenus, tspecies, tsciName, tauthority = "","","","","","","","",""
		if irmng_data != {}:
			if 'kingdom' in irmng_data:
				tkingdom = irmng_data['kingdom']
			if 'phylum' in irmng_data:
				tphylum = irmng_data['phylum']
			if 'class' in irmng_data:
				tclass = irmng_data['class']
			if 'order' in irmng_data:
				torder = irmng_data['order']
			if 'family' in irmng_data:
				tfamily = irmng_data['family']
			if 'genus' in irmng_data:
				tgenus = irmng_data['genus']
			if irmng_data["rank"] == "Species":
				if 'scientificname' in irmng_data:
					tspecies = irmng_data['scientificname'].split()[1]
					tsciName=irmng_data['scientificname']
				if 'authority' in irmng_data:
					tauthority = irmng_data['authority']
		
		out_dict={
			"Kingdom": tkingdom,
			"Phylum": tphylum,
			"Class": tclass,
			"taxOrder": torder,
			"Family": tfamily,
			"Genus": tgenus,
			"Species": tspecies,
			"ScientificName": tsciName,
			"Authority": tauthority
			}
		return out_dict
	
	# function for retrieving the taxonomic path from PESI
	def pesiTaxpath(self):
		"""
		Function for getting the taxnonomic data from the PESI response.
		Returns a dictionary of strings.
		"""

		pesi_data = self.response_dict["PESI"]
		tkingdom, tphylum, tclass, torder, tfamily, tgenus, tspecies, tsciName, tauthority = "","","","","","","","",""
		if pesi_data != {}:
			if 'kingdom' in pesi_data:
				tkingdom = pesi_data['kingdom']
			if 'phylum' in pesi_data:
				tphylum = pesi_data['phylum']
			if 'class' in pesi_data:
				tclass = pesi_data['class']
			if 'order' in pesi_data:
				torder = pesi_data['order']
			if 'family' in pesi_data:
				tfamily=pesi_data['family']
			if 'genus' in pesi_data:
				tgenus = pesi_data['genus']
			if pesi_data["rank"] == "Species":
				if 'scientificname' in pesi_data:
					tspecies = pesi_data['scientificname'].split()[1]
					tsciName = pesi_data['scientificname']
				if 'authority' in pesi_data:
					tauthority = pesi_data['authority']
		
		out_dict={
			"Kingdom": tkingdom,
			"Phylum": tphylum,
			"Class": tclass,
			"taxOrder": torder,
			"Family": tfamily,
			"Genus": tgenus,
			"Species": tspecies,
			"ScientificName": tsciName,
			"Authority": tauthority
			}
		return out_dict

	def wormsTaxpath(self):
		"""
		Function for getting the taxnonomic data from the WORMS response.
		Returns a dictionary of strings.
		"""

		worms_data = self.response_dict["WORMS"]
		tkingdom, tphylum, tclass, torder, tfamily, tgenus, tspecies, tsciName, tauthority = "","","","","","","","",""
		if worms_data != {}:
			if 'kingdom' in worms_data:
				tkingdom = worms_data['kingdom']
			if 'phylum' in worms_data:
				tphylum = worms_data['phylum']
			if 'class' in worms_data:
				tclass = worms_data['class']
			if 'order' in worms_data:
				torder = worms_data['order']
			if 'family' in worms_data:
				tfamily = worms_data['family']
			if 'genus' in worms_data:
				tgenus = worms_data['genus']
			if worms_data["rank"] == "Species":
				if 'scientificname' in worms_data:
					tspecies = worms_data['scientificname'].split()[1]
					tsciName = worms_data['scientificname']
				if 'authority' in worms_data:
					tauthority = worms_data['authority']
		
		out_dict={
			"Kingdom": tkingdom,
			"Phylum": tphylum,
			"Class": tclass,
			"taxOrder": torder,
			"Family": tfamily,
			"Genus": tgenus,
			"Species": tspecies,
			"ScientificName": tsciName,
			"Authority": tauthority
			}
		return out_dict

	def irmngTraits(self):
		"""
		Function for getting the trait data from the IRMNG response.
		Returns a dictionary of booleans.
		"""

		irmng_data = self.response_dict["IRMNG"]

		isMarine, isBrackish, isFreshwater, isTerrestrial, isAllWater, isMarineFresh, isExtinct = None,None,None,None,None,None,None
		if 'isMarine' in irmng_data:
			isMarine = irmng_data['isMarine']
		if 'isBrackish' in irmng_data:
			isBrackish = irmng_data['isBrackish']
		if 'isFreshwater' in irmng_data:
			isFreshwater = irmng_data['isFreshwater']
		if 'isTerrestrial' in irmng_data:
			isTerrestrial = irmng_data['isTerrestrial']
		if 'isExtinct' in irmng_data:
			isExtinct = irmng_data['isExtinct']

		# if the species lives in all water habitats
		if isMarine == 1 and isFreshwater == 1 and isBrackish == 1:
			isAllWater = 1
		else:
			isAllWater = 0
		# if the species lives in marine and freshwater habitats
		if isMarine == 1 and isFreshwater == 1:
			isMarineFresh = 1
		else:
			isMarineFresh = 0
		
		out_dict={
			"isMarine": isMarine,
			"isBrackish": isBrackish,
			"isFreshwater": isFreshwater,
			"isTerrestrial": isTerrestrial,
			"isAllWater": isAllWater,
			"isMarineFresh": isMarineFresh,
			"isExtinct": isExtinct
			}
		return out_dict

	def wormsTraits(self):
		"""
		Function for getting the trait data from the WORMS response.
		Returns a dictionary of booleans.
		"""

		worms_data = self.response_dict["WORMS"]

		isMarine, isBrackish, isFreshwater, isTerrestrial, isAllWater, isMarineFresh, isExtinct = None,None,None,None,None,None,None
		if 'isMarine' in worms_data:
			isMarine = worms_data['isMarine']
		if 'isBrackish' in worms_data:
			isBrackish = worms_data['isBrackish']
		if 'isFreshwater' in worms_data:
			isFreshwater = worms_data['isFreshwater']
		if 'isTerrestrial' in worms_data:
			isTerrestrial = worms_data['isTerrestrial']
		if 'isExtinct' in worms_data:
			isExtinct = worms_data['isExtinct']
		
		# if the species lives in all water habitats
		if isMarine == 1 and isFreshwater == 1 and isBrackish == 1:
			isAllWater = 1
		else:
			isAllWater = 0
		# if the species lives in marine and freshwater habitats
		if isMarine == 1 and isFreshwater == 1:
			isMarineFresh = 1
		else:
			isMarineFresh = 0

		out_dict={
			"isMarine": isMarine,
			"isBrackish": isBrackish,
			"isFreshwater": isFreshwater,
			"isTerrestrial": isTerrestrial,
			"isAllWater": isAllWater,
			"isMarineFresh": isMarineFresh,
			"isExtinct": isExtinct
			}
		return out_dict

	def allVernaculars(self):
		ver_dict = {
			"Wikipedia": self.response_dict["Wikipedia-Vernaculars"],
			"GBIF": self.response_dict["GBIF-Vernaculars"],
			"WORMS": self.response_dict["WORMS-Vernaculars"],
		}
		return ver_dict

	def combineDictionaries(self):
		"""
		Function for combining all relevant dictionaries into one single dict.
		Returns a dictionary of dictionaries.
		"""
		combined_dict = {
			"Taxonomy": {
				"WORMS": self.wormsTaxpath(),
				"IRMNG": self.irmngTaxpath(),
				"GBIF": self.gbifTaxpath(),
				"PESI": self.pesiTaxpath(),
				"Vernaculars": self.allVernaculars()
			},
			"Traits": {
				"WORMS": self.wormsTraits(),
				"IRMNG": self.irmngTraits()
			},
			"IDs": self.id_dict
		}
		return(combined_dict)

class ResolveData:
	def __init__(self, sciName: str, accessionNumber: str = ""):
		search = GetInformation(sciName, accessionNumber)
		self.data_dict = search.combineDictionaries()
	
	def resolveTaxonomy(self):
		"""
		Resolves the taxonomy dictionary by iterating over the keys one by one.
		If no WORMS data is available, the next key is tried, and so forth.
		Returns a dictionary with the resolved data. Dict is empty if no data is available.
		"""
		tax_dict = self.data_dict["Taxonomy"]
		# get keys from the first subdictionary (excluding 'Vernaculars')
		keys = [key for key in tax_dict['WORMS'].keys() if key != 'Vernaculars']
		main_keys = [key for key in tax_dict.keys() if key != 'Vernaculars']

		# initialize output dictionary
		out_dict = {}
		# iterate through taxonomy dict
		for main_key in main_keys:
			# set output dict values to values of current sub dictionary
			out_dict = {key: tax_dict[main_key][key] for key in keys}
			# if all values of current dictionary are empty, continue. else, break out of for loop
			if all(out_dict.values()):
				break

		return(out_dict)
	
	def resolveVernaculars(self):
		ver_dict = self.data_dict["Taxonomy"]["Vernaculars"]
		
		out_dict = {"English": "", "German": ""}
		for key in ver_dict:
			if ver_dict[key]["English"] != '' and out_dict["English"] == '':
				out_dict["English"] = ver_dict[key]["English"]
			if ver_dict[key]["German"] != '' and out_dict["German"] == '':
				out_dict["German"] = ver_dict[key]["German"]
		
		return out_dict


	def resolveTraits(self):
		"""
		Resolves the trait dictionary by comparing all subdictionaries and assigning values according to the following rules:
		If all values for a trait are the same, set that value (0 or 1)
		If all values for a trait are None, set 2
		If conflicts for a trait can not be resolved, set 3
		Returns a dictionary with the resolved values.
		"""

		# extract trait directory
		trait_dict = self.data_dict["Traits"]
		# get keys from the first subdictionary (assuming all have the same keys)
		keys = list(next(iter(trait_dict.values())).keys())
		
		# initialize output dictionary
		out_dict = {}
		 # iterate over each key
		for key in keys:
			# get values for the current key from each subdictionary
			values = [subdict.get(key) for subdict in trait_dict.values()]
			# filter out all None values
			filtered_values = [value for value in values if value is not None]

			# check conditions and set output value accordingly
			if not filtered_values:
				# set to 2 if all values are None
				out_dict[key] = 2
			elif len(set(filtered_values)) == 1:
				# set to 1 or 0 if all values are the same
				out_dict[key] = filtered_values[0]
			else:
				# set to 3 in any other case
				out_dict[key] = 3
		
		# if the combined habitat values are not matching, set them again
		if out_dict["isAllWater"] == 3:
			if out_dict["isMarine"] == 1 and out_dict["isFreshwater"] == 1 and out_dict["isBrackish"] == 1:
				out_dict["isAllWater"] = 1
			else:
				out_dict["isAllWater"] = 0
		if out_dict["isMarineFresh"] == 3:
			if out_dict["isMarine"] == 1 and out_dict["isFreshwater"] == 1:
				out_dict["isMarineFresh"] = 1
			else:
				out_dict["isMarineFresh"] = 0

		return out_dict
	
	def combineDictionaries(self):
		resolved_dict = {
			"Taxonomy": self.resolveTaxonomy(),
			"Vernaculars": self.resolveVernaculars(),
			"Traits": self.resolveTraits(),
			"IDs": self.data_dict["IDs"]
		}
		return(resolved_dict)

if __name__=="__main__":
	# cannot deal with subspecies right now
	testSciName = "Oncorhynchus mykiss"
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
	print(f"Test Species: {testValues['sciNameList'][0]}")

	getScientificNames(testValues["accessionNumberList"])

	#data = ResolveData(testSciName)
	#print(data.combineDictionaries())
	#print("\n")

	#search = GetInformation(testSciName)
	#print(search.combineDictionaries())
	#print("\n")









