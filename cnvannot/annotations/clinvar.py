import csv
from intervaltree import IntervalTree

from cnvannot.common.serialization import *

# Clinvar
def clinvar_load():
    chr_dict = {}
    base_file = 'clinvar.tsv'
    base_path = os.path.join(Common.data_path, base_file)

    if serialization_is_serialized(base_file):
        return serialization_deserialize(base_file)

    with open(base_path) as csvfile:
        csvreader = csv.reader(csvfile,delimiter="\t")
        for row in csvreader:
                
            #Interval tree object
            chrom = "chr" + row[0]
            start =  int(row[1])
            stop = int(row[2])
            var_type="loss"
            
            if row[3] == "1":
                var_type="gain"
			
            itobj = {'chr':chrom, 'start': start, 'stop': stop, 'var_type': var_type, 'gene': row[4], 'isPatho':row[5], 'phenotype':row[6], 'dbs':row[7]  }
                
            #Add to interval trees, 1 interval tree by chromosome
            if chrom not in chr_dict:
                chr_dict[chrom] = IntervalTree()
            try:
                chr_dict[chrom][start:stop] = itobj
            except ValueError:
                pass

    serialization_serialize(chr_dict, base_file)
    return chr_dict
