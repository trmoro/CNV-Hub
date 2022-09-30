import os
import csv
import json

from cnvannot.annotations.dgv import dgv_gold_load
from cnvannot.annotations.encode import encode_load
from cnvannot.annotations.omim import omim_morbid_genes_load
from cnvannot.annotations.refseq import refseq_load
from cnvannot.annotations.special.france_incomplete_penetrance import france_incomplete_penetrance_load
from cnvannot.annotations.ucsc import ucsc_get_annotation_link
from cnvannot.annotations.decipher import decipher_load

from cnvannot.common.coordinates import coordinates_from_string

from cnvannot.queries.basic_queries import *
from cnvannot.queries.interpretation import interpretation_get
from cnvannot.queries.specific_queries import dgv_gold_overlap_count_1_percent, exc_overlaps_70_percent, omim_match_organ

from cloud.aws import getS3, readFile, setFile

dgv_db = dgv_gold_load()
refseq_db = refseq_load()
omim_mg_db = omim_morbid_genes_load()
encode_db = encode_load()
france_inc_pen_db = france_incomplete_penetrance_load()
decipher_db = decipher_load()

#XCNV MVP Score Interpretation
def xcnv_interpretation_from_score(score: float) -> str:
    if score < 0.25:
        return 'benign'
    elif score < 0.5:
        return 'likely benign'
    elif score < 0.75:
        return 'likely pathogenic'
    else:
        return 'pathogenic'
        
#Read XCNV
def read_xcnv_scores():
    print("XCNV Score : ")
    scores = []
    with open("data.output.csv") as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            scores.append(float(row["MVP_score"]))
    print(scores)
    return scores

#Search
def search(str_query: str, organ: str, score: float):
	
    print("Searching : " + str_query)
	
    query = coordinates_from_string(str_query)
    ucsc_url = str(ucsc_get_annotation_link(query)['ucsc']['link'])
    cnv_len = query.end - query.start
    cnv_type = str.upper(query.type)
    exclude_overlaps = exc_overlaps_70_percent(encode_db, query)
    
    gene_overlap_count = query_overlap_count(refseq_db, query)
    morbid_gene_overlap_count = query_overlap_count(omim_mg_db, query)
    decipher_overlap_count = 0
    
    gene_overlaps = []
    morbid_gene_overlaps = []
    decipher_overlaps = []

    for g in query_overlaps_get(refseq_db, query):
        gene_overlaps.append(g[2])
	
    for m in query_overlaps_get(omim_mg_db, query):
        morbid_gene_overlaps.append(m[2])
        
    for d in query_overlaps_get(decipher_db, query):
        if d.data['var_type'] == query.type:
            decipher_overlaps.append(d[2])
    decipher_overlap_count = len(decipher_overlaps)

    organ_match_count = omim_match_organ(omim_mg_db, query, organ)
    dgv_gold_cnv_overlap_count = dgv_gold_overlap_count_1_percent(dgv_db, query)

    xcnv_res = score
    xcnv_res_interpretation = xcnv_interpretation_from_score(xcnv_res)

    synth_interpretation = 'Interpretation suggestion(s): ' + interpretation_get(query,
                                                                                 xcnv_res,
                                                                                 exclude_overlaps,
                                                                                 gene_overlap_count,
                                                                                 morbid_gene_overlap_count,
                                                                                 dgv_gold_cnv_overlap_count,
                                                                                 query.type,
                                                                                 france_inc_pen_db,
                                                                                 organ_match_count
                                                                                 )[0:-68]

    title = str_query
    if organ != "":
	    title += " (" + organ + ")"
    return {		'title': title,
					'ucsc_url': ucsc_url,
                    'cnv_len': cnv_len,
                    'cnv_type': cnv_type,
                    'exc_overlaps': exclude_overlaps,
                    'gene_overlap_count': gene_overlap_count,
                    'morbid_gene_overlap_count': morbid_gene_overlap_count,
                    'morbid_gene_pheno_overlap_count': organ_match_count,
                    'decipher_overlap_count': decipher_overlap_count,
                    'dgv_gold_cnv_overlap_count': dgv_gold_cnv_overlap_count,
                    'xcnv_res': xcnv_res,
                    'xcnv_res_interpretation': xcnv_res_interpretation,
                    'synth_interpretation': synth_interpretation,
                    'gene_overlaps' : gene_overlaps,
                    'morbid_gene_overlaps': morbid_gene_overlaps,
                    'decipher_overlaps': decipher_overlaps}
                    
#Test
def test_xcnv():
    os.system("./XCNV.R XCNV/example_data/1.bed")
    with open("XCNV/example_data/1.output.csv") as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            print(row)

#Download of S3
def downloadS3(s3Path,localPath):
    conn, bucket = getS3()
    data = readFile(bucket,s3Path).decode("utf-8") 
    filedata = open(localPath,'w',errors='ignore')
    filedata.write(str(data) )
    filedata.close()
    return data

#Save to S3
def saveS3(localPath,s3Path):
    conn, bucket = getS3()
    filedata = open(localPath,'r',errors='ignore')
    setFile(bucket,s3Path,filedata.read() )
    filedata.close()

#Process
def process(filepath):
    userID = filepath.split("/")[0]
    downloadS3(filepath,"data.bed")
    ref = "hg19"
    
    #XCNV
    os.system("./XCNV.R data.bed")
    saveS3("data.output.csv",filepath + ".output.xcnv")
    scores = read_xcnv_scores()
    
    #CNV Hub
    res = []
    f = open("data.bed","r")
    i = 0
    for l in f:
        d = l.split("\t")
        query = ref + ":chr" + d[0] + ":" + d[1] + "-" + d[2] + ":" + d[3]
        res.append(search(query,"",scores[i]))
        i+=1
    f.close()
    print("CNV-Hub Result : ")
    print(res)
    print("--------------------")
    conn, bucket = getS3()
    setFile(bucket,filepath + ".output.json",json.dumps(res))
