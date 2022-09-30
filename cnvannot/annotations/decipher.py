import csv
from intervaltree import IntervalTree

from cnvannot.common.serialization import *

# Decipher
def decipher_load():
    chr_dict = {}
    base_file = 'decipher.csv'
    base_path = os.path.join(Common.data_path, base_file)

    if serialization_is_serialized(base_file):
        return serialization_deserialize(base_file)

    with open(base_path) as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
			
            #Get data from each row
            data = row[2].split(" ")
			
            #Only gain and losses
            if data[2] == "Duplication" or data[2] == "Deletion":				
                var_type = "loss"
                if data[2] == "Duplication":
                    var_type = "gain"
                chrom = "chr" + data[0].split(":")[0]
                start = int(data[0].split(":")[1].split("-")[0])
                stop = int(data[0].split(":")[1].split("-")[1])
                heterozygous = row[4]
                pathologic = row[5]
                description = row[6]
                
                #Interval tree object
                itobj = {'chr':chrom, 'start': start, 'stop': stop, 'var_type': var_type, 'heterozygous':heterozygous, 'pathologic':pathologic, 'description':description  }
                
                #Add to interval trees, 1 interval tree by chromosome
                if chrom not in chr_dict:
                    chr_dict[chrom] = IntervalTree()
                try:
                    chr_dict[chrom][start:stop] = itobj
                except ValueError:
                    pass

    serialization_serialize(chr_dict, base_file)
    return chr_dict
