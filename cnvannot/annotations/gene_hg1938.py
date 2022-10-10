from cnvannot.common.serialization import *
from intervaltree import IntervalTree

# Load
def gene_hg1938_load():
	chr_dict = {}

	hg19_path = os.path.join(Common.data_path, 'hg19.tsv')
	hg38_path = os.path.join(Common.data_path, 'hg38.tsv')
	hg1918 = 'gene_hg1938.bin'

	if serialization_is_serialized(hg1918):
		return serialization_deserialize(hg1918)

	with open(hg19_path) as f:
		lines = f.readlines()
		for line in lines:
			parts = line.split('\t')
			chrom = parts[2]
			start = int(parts[4])
			stop = int(parts[5])
			gene = parts[12]
			
			#Add to interval trees, 1 interval tree by chromosome
			itobj = {'chr':chrom, 'start': start, 'stop': stop, 'gene': gene, 'ref': 'hg19' }
			if chrom not in chr_dict:
				chr_dict[chrom] = IntervalTree()
			try:
				chr_dict[chrom][start:stop] = itobj
			except ValueError:
				pass
	
	with open(hg38_path) as f:
		lines = f.readlines()
		for line in lines:
			parts = line.split('\t')
			chrom = parts[2]
			start = int(parts[4])
			stop = int(parts[5])
			gene = parts[12]
			
			#Add to interval trees, 1 interval tree by chromosome
			itobj = {'chr':chrom, 'start': start, 'stop': stop, 'gene': gene, 'ref': 'hg38' }
			if chrom not in chr_dict:
				chr_dict[chrom] = IntervalTree()
			try:
				chr_dict[chrom][start:stop] = itobj
			except ValueError:
				pass		

	serialization_serialize(chr_dict, hg1918)

	return chr_dict
