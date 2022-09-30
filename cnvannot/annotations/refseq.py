from intervaltree import IntervalTree
from cnvannot.common.serialization import *


# RefSeq DB

def refseq_load():
    chr_dict = {}

    refseq_base_file = 'refGene.txt'
    refseq_base_path = os.path.join(Common.data_path, refseq_base_file)

    if serialization_is_serialized(refseq_base_file):
        return serialization_deserialize(refseq_base_file)

    with open(refseq_base_path) as f:
        for line in f:
            parts = line.split('\t')
            chrom = parts[2]
            start = int(parts[4])
            stop = int(parts[5])
            name1 = parts[1].lower()
            name2 = parts[12].lower()

            if chrom not in chr_dict:
                # Add new interval tree as value
                chr_dict[chrom] = IntervalTree()
            pass

            try:
                chr_dict[chrom][start:stop] = {'chr': chrom, 'start': start, 'stop': stop,
                                               'name1': name1, 'name2': name2}
            except ValueError:
                pass

    serialization_serialize(chr_dict, refseq_base_file)

    return chr_dict
