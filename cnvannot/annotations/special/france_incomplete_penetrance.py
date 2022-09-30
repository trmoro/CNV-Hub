from intervaltree import IntervalTree

from cnvannot.common.serialization import *


# FRANCE INCOMPLETE PENETRANCE (PARTIAL)

def france_incomplete_penetrance_load():
    chr_dict = {}

    fip_base_file = 'france_incomplete_penetrance.hg19.txt'
    fip_base_path = os.path.join(Common.data_path, fip_base_file)

    if serialization_is_serialized(fip_base_file):
        return serialization_deserialize(fip_base_file)

    with open(fip_base_path) as f:
        for line in f:
            if line.strip() == '':
                continue
            parts = line.split('\t')
            chrom = parts[0]
            start = int(parts[1])
            stop = int(parts[2])
            cnv_type = parts[3]
            desc = parts[4]

            if chrom not in chr_dict:
                # Add new interval tree as value
                chr_dict[chrom] = IntervalTree()

            try:
                chr_dict[chrom][start:stop] = {'chr': chrom, 'start': start, 'stop': stop,
                                               'var_type': cnv_type, 'desc': desc}
            except ValueError:
                pass

    serialization_serialize(chr_dict, fip_base_file)

    return chr_dict
