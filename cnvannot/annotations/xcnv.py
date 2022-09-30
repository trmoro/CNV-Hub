from time import sleep

import docker
from docker.errors import DockerException
from intervaltree import IntervalTree

# XCNV DB (ONLY hg19)
from cnvannot.annotations.base import base_coordinates_annotation


def xcnv_is_avail() -> bool:
    try:
        client = docker.from_env()
        images = client.images
    except DockerException:
        return False

    for img in images.list():
        if "xcnvcl" in str(img):
            return True

    return False


def xcnv_predict(queries) -> list:
    xcnv_query = ''
    for q in queries:
        if q.ref != 'hg19':
            raise Exception('XCNV works only on hg19')
        xcnv_query += q.chr + ':' + str(q.start) + '-' + str(q.end) + ':' + q.type + ','
    xcnv_query = xcnv_query[:-1]

    client = docker.from_env()

    cont = client.containers.run('xcnvcl', xcnv_query, detach=True)

    timeout = True

    chr_dict = {}

    # Wait loop.
    for i in range(300):
        curr_out = cont.logs()
        if len(curr_out) > 0:
            # Output received.
            timeout = False
            possible_prediction = curr_out.split()[-1]
            prediction_parts = possible_prediction.decode('ascii').split(',')
            sanity_check_count = int(prediction_parts[0])
            if len(prediction_parts) != (sanity_check_count + 2):
                raise Exception("X-CNV sanity check failed!")

            if prediction_parts[-1].strip() != 'EOF':
                raise Exception("X-CNV sanity check failed!")

            for p in prediction_parts[1:-1]:
                p_parts = p.split('_')
                chrom = 'chr' + p_parts[0]
                if chrom not in chr_dict:
                    chr_dict[chrom] = IntervalTree()
                chr_dict[chrom][int(p_parts[1]):int(p_parts[2])] = float(p_parts[-1])

            break

        sleep(0.5)

    if timeout:
        raise Exception("Timeout Exception")

    ret_dict_list = []
    for query in queries:
        ret_dict = base_coordinates_annotation(query)
        itree: IntervalTree = chr_dict[query.chr]
        for interval in itree[query.start:query.end]:
            ret_dict["xcnv"] = {'prediction': interval.data}
        ret_dict_list.append(ret_dict)

    return ret_dict_list


def xcnv_interpretation_from_score(score: float) -> str:
    if score < 0.25:
        return 'benign'
    elif score < 0.5:
        return 'likely benign'
    elif score < 0.75:
        return 'likely pathogenic'
    else:
        return 'pathogenic'
