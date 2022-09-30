from intervaltree import IntervalTree

from cnvannot.common.serialization import *


# DGV DB

def dgv_gold_load():
    chr_dict = {}

    dgv_base_file = 'DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3'
    dgv_base_path = os.path.join(Common.data_path, dgv_base_file)

    if serialization_is_serialized(dgv_base_file):
        return serialization_deserialize(dgv_base_file)

    with open(dgv_base_path) as f:
        for line in f:
            parts = line.split('\t')
            chrom = parts[0]
            start = int(parts[3])  # TODO: A Verifier
            stop = int(parts[4])  # TODO: A VERIFIER
            dat = parts[8]
            dat_parts = dat.split(';')
            var_type = ''
            freq_percent = ''
            for dp in dat_parts:
                if dp.startswith('variant_sub_type'):
                    var_type = dp.split('=')[1]
                if dp.startswith('Frequency'):
                    freq_percent = float(dp.split('=')[1][:-1])
            if var_type == 'Gain':
                var_type = 'GAIN'
            elif var_type == 'Loss':
                var_type = 'LOSS'
            if var_type == '':
                print("Cannot parse line")

            if chrom not in chr_dict:
                # Add new interval tree as value
                chr_dict[chrom] = IntervalTree()

            try:
                chr_dict[chrom][start:stop] = {'chr': chrom, 'start': start, 'stop': stop,
                                               'var_type': var_type, 'freq': freq_percent}
            except ValueError:
                pass

    serialization_serialize(chr_dict, dgv_base_file)

    return chr_dict
