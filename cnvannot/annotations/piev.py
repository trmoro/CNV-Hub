import csv
from intervaltree import IntervalTree

from cnvannot.common.serialization import *

# Piev
def piev_load():
    chr_dict = {}
    base_file = 'piev.csv'
    base_path = os.path.join(Common.data_path, base_file)

    if serialization_is_serialized(base_file):
        return serialization_deserialize(base_file)

    with open(base_path) as csvfile:
        csvreader = csv.reader(csvfile,delimiter=";")
        for row in csvreader:
                
            #Interval tree object
            chrom = row[0]
            start =  int(row[1])
            stop = int(row[2])
            var_type=row[3]
			
            itobj = {'chr':chrom, 'start': start, 'stop': stop, 'var_type': var_type}
                
            #Add to interval trees, 1 interval tree by chromosome
            if chrom not in chr_dict:
                chr_dict[chrom] = IntervalTree()
            try:
                chr_dict[chrom][start:stop] = itobj
            except ValueError:
                pass

    serialization_serialize(chr_dict, base_file)
    return chr_dict
