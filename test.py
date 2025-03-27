from Bio import Entrez
import xml.etree.ElementTree as ET

def get_scientific_names(accession_numbers):
    Entrez.email = "your_email@example.com"
    scientific_names = {}
    
    for accession in accession_numbers:
        try:
            # Try using esummary first
            handle = Entrez.esummary(db="assembly", id=accession)
            record = Entrez.read(handle)
            if record[0].get("Organism_Name"):
                scientific_name = record[0]["Organism_Name"]
            else:
                # If esummary doesn't work, try parsing XML manually
                handle = Entrez.efetch(db="assembly", id=accession, rettype="gb", retmode="xml")
                xml_content = handle.read()
                root = ET.fromstring(xml_content)
                # You'll need to navigate the XML tree here to find the organism name
                # This part depends on the structure of the XML response
                scientific_name = "Not Found"
            
            scientific_names[accession] = scientific_name
        except Exception as e:
            print(f"Error processing accession {accession}: {e}")
    
    return scientific_names

# Example usage
accession_numbers = ["GCF_020740795.1", "GCA_033628465.1"]
scientific_names_dict = get_scientific_names(accession_numbers)

for accession, name in scientific_names_dict.items():
    print(f"{accession}: {name}")
