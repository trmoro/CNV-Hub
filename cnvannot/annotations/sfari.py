from cnvannot.common.serialization import *


# Sfari
def sfari_load():
	gene_dict = {}

	sfari_base_file = 'sfari_simplified.csv'
	sfari_base_path = os.path.join(Common.data_path, sfari_base_file)

	if serialization_is_serialized(sfari_base_file):
		return serialization_deserialize(sfari_base_file)

	with open(sfari_base_path) as f:
		lines = f.readlines()
		for line in lines:
			parts = line.split(',')
			gene = parts[0]
			gene_score = parts[1]
			syndromic = parts[2]
			eagle = parts[3]
			reports = parts[4]
			
			if gene_score != "":
				gene_score = int(gene_score)
			else:
				gene_score = -1
				
			if syndromic != "":
				syndromic = int(syndromic)
			else:
				syndromic = -1
				
			if eagle != "":
				eagle = float(eagle)
			else:
				eagle = -1
				
			if reports != "":
				reports = int(reports)
			else:
				reports = -1
				
			gene_dict[gene] = { "gene_score" : gene_score, "syndromic" : syndromic, "eagle" : eagle, "reports": reports}
				

	serialization_serialize(gene_dict, sfari_base_file)

	return gene_dict
