from intervaltree import IntervalTree
from cnvannot.common.serialization import *


# Encode DB

def encode_load():
    chr_dict = {}

    encode_base_file = 'consensusBlacklist.bed'
    encode_base_path = os.path.join(Common.data_path, encode_base_file)

    if serialization_is_serialized(encode_base_file):
        return serialization_deserialize(encode_base_file)

    with open(encode_base_path) as f:
        for line in f:
            parts = line.split('\t')
            chrom = parts[0]
            start = int(parts[1])
            stop = int(parts[2])
            dat = parts[3]

            if chrom not in chr_dict:
                # Add new interval tree as value
                chr_dict[chrom] = IntervalTree()

            try:
                chr_dict[chrom][start:stop] = {'exclusion_reason': dat}
            except ValueError:
                pass

    serialization_serialize(chr_dict, encode_base_file)

    return chr_dict
