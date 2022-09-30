from cnvannot.common.serialization import *

# Cast Float
def castFloat(value : str) -> float:
	if value != "-" and value != "":
		return float(value)
	else:
		return -1


# Load
def gene_stat_load():
	gene_dict = {}

	stat_base_file = 'gene_stat.csv'
	stat_base_path = os.path.join(Common.data_path, stat_base_file)

	if serialization_is_serialized(stat_base_file):
		return serialization_deserialize(stat_base_file)

	with open(stat_base_path) as f:
		lines = f.readlines()
		for line in lines:
			parts = line.replace("\n","").split(';')
			gene = parts[0]
			pLI = castFloat(parts[1])
			LOEUF = castFloat(parts[2])
			sHet = castFloat(parts[3])
			pHaplo = castFloat(parts[4])
			pTriplo = castFloat(parts[5])
				
			gene_dict[gene] = { "pLI" : pLI, "LOEUF" : LOEUF, "sHet" : sHet, "pHaplo": pHaplo, "pTriplo" : pTriplo}
				

	serialization_serialize(gene_dict, stat_base_file)

	return gene_dict
