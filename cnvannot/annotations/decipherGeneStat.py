import csv
from intervaltree import IntervalTree

from cnvannot.common.serialization import *

# Cast Float
def castFloat(value : str) -> float:
    if value != "-":
	    return float(value)
    else:
	    return -1

# Decipher Gene Stat
def decipherGeneStat_load():
    chr_dict = {}
    base_file = 'decipherGeneStat.csv'
    base_path = os.path.join(Common.data_path, base_file)

    if serialization_is_serialized(base_file):
        return serialization_deserialize(base_file)

    with open(base_path) as csvfile:
        csvreader = csv.reader(csvfile, delimiter=';')
        for row in csvreader:
            #print(row)
            chrom = "chr" + row[1].split(":")[0]
            start = int(row[1].split(":")[1].split("-")[0])
            stop = int(row[1].split(":")[1].split("-")[1])
            gene = row[0]
            pLI = castFloat(row[2])
            LOEUF = castFloat(row[3])
            sHet = castFloat(row[4])
            pHaplo = castFloat(row[5])
            pTriplo = castFloat(row[6])
                
            #Interval tree object
            itobj = {'chr':chrom, 'start': start, 'stop': stop, 'gene' : gene, 'pLI' : pLI, 'LOEUF': LOEUF, 'sHet': sHet, 'pHaplo': pHaplo, 'pTriplo': pTriplo }
            #print(itobj)
                
            #Add to interval trees, 1 interval tree by chromosome
            if chrom not in chr_dict:
                chr_dict[chrom] = IntervalTree()
            try:
                chr_dict[chrom][start:stop] = itobj
            except ValueError:
                pass

    serialization_serialize(chr_dict, base_file)
    return chr_dict
