import os
import csv
import json

from cnvannot.annotations.dgv import dgv_gold_load
from cnvannot.annotations.encode import encode_load
from cnvannot.annotations.omim import omim_morbid_genes_load
from cnvannot.annotations.refseq import refseq_load
from cnvannot.annotations.ucsc import ucsc_get_annotation_link
from cnvannot.annotations.decipher import decipher_load
from cnvannot.annotations.clinvar import clinvar_load
from cnvannot.annotations.sfari import sfari_load
from cnvannot.annotations.gene_stat import gene_stat_load
from cnvannot.annotations.dijon import dijon_load
from cnvannot.annotations.gene_hg1938 import gene_hg1938_load

from cnvannot.common.coordinates import coordinates_from_string

from cnvannot.queries.basic_queries import query_overlaps_get, compute_overlaps, filter_hg1938
from cnvannot.queries.specific_queries import exc_overlaps_70_percent, omim_get_organs, inject_sfari, inject_gene_stat
from cnvannot.queries.scoring import compute_scores

from cloud.aws import getS3, setFile, deleteFile, readFile

#Databases
dgv_db = dgv_gold_load()
refseq_db = refseq_load()
omim_mg_db = omim_morbid_genes_load()
encode_db = encode_load()
decipher_db = decipher_load()
clinvar_db = clinvar_load()
sfari_db = sfari_load()
stat_db = gene_stat_load()
dijon_db = dijon_load()
hg1938_db = gene_hg1938_load()	

#Save result
def save(data : str, query: str, organ : str):
	filename = "cnv/" + query
	if organ != "":
		filename += "_" + organ
	conn, bucket = getS3()
	setFile(bucket,filename + ".json",json.dumps(data))
	conn.close()

#XCNV
def xcnv(str_query: str):
	
	print("Processing XCNV - " + str_query)
		
	#Process query to .bed
	data = str_query.split(":")
	str_data = data[1].replace("chr","") + "\t" + data[2].replace("-","\t") + "\t" + data[3]
	
	#Write file
	f = open("/tmp/tmp.bed","w")
	f.write(str_data)
	f.close()
		
	#Execute XCNV and convert result to regular dict
	os.system("./XCNV/bin/XCNV /tmp/tmp.bed")
	with open("/tmp/tmp.output.csv") as csvfile:
		csvreader = csv.DictReader(csvfile)
		for row in csvreader:
			return json.loads(json.dumps(row))

#Generate CNV Object
def gen_cnv(title,ucsc_url,scores,cnv_len,query,cnv_type,exclude_overlaps,morbid_gene_overlaps,decipher_overlaps,clinvar_overlaps,dgv_overlaps,dijon_overlaps,ref_overlaps,xcnv_data):
	return {'title': title,
	'ucsc_url': ucsc_url,
	'score_data': scores,
	'cnv_len': cnv_len,
	'cnv_ref': query.ref,
	'cnv_chr': query.chr,
	'cnv_start': query.start,
	'cnv_end': query.end,
	'cnv_type': cnv_type,
	'exc_overlaps': exclude_overlaps,
	'morbid_gene_overlaps': morbid_gene_overlaps,
	'decipher_overlaps': decipher_overlaps,
	'clinvar_overlaps': clinvar_overlaps,
	'dgv_overlaps': dgv_overlaps,
	'dijon_overlaps': dijon_overlaps,
	'gene_overlaps': ref_overlaps,
	'xcnv': xcnv_data}

#Search
def search(str_query: str, organ: str):
	
    print("Searching : " + str_query)
	
    query = coordinates_from_string(str_query)
    ucsc_url = str(ucsc_get_annotation_link(query)['ucsc']['link'])
    cnv_len = query.end - query.start
    cnv_type = str.upper(query.type)
    exclude_overlaps = exc_overlaps_70_percent(encode_db, query)
    
    #Compute overlaps on ref database
    ref_overlaps = compute_overlaps(hg1938_db, query)
    ref_overlaps = filter_hg1938(ref_overlaps, query.ref, ref_overlaps)
    	
	#Compute Overlaps on OMIM, Decipher, Decipher Stat, Clinvar, DGV and CHU DIJON, with overlaps rate
    morbid_gene_overlaps = filter_hg1938(compute_overlaps(omim_mg_db, query), query.ref, ref_overlaps)
    decipher_overlaps = compute_overlaps(decipher_db, query)
    clinvar_overlaps = filter_hg1938(compute_overlaps(clinvar_db, query), query.ref, ref_overlaps)
    dgv_overlaps = compute_overlaps(dgv_db,query,False,0.1)   
    dijon_overlaps = compute_overlaps(dijon_db, query)

    #Inject Sfari and Gene stat in OMIM
    inject_sfari(morbid_gene_overlaps,sfari_db)
    inject_gene_stat(morbid_gene_overlaps,stat_db)
    
    #Inject Sfari and Gene stat in RefGene
    inject_sfari(ref_overlaps,sfari_db)
    inject_gene_stat(ref_overlaps,stat_db)
     
    #Organ list
    organ_list = omim_get_organs(omim_mg_db, query)
           
	#Compute Scores
    scores = compute_scores(cnv_len,dgv_overlaps,clinvar_overlaps,decipher_overlaps,dijon_overlaps)
    #print(scores)
	
	#Title
    title = str_query
    if organ != "":
	    title += " (" + organ + ")"
	    
	#Save without XCNV
    print("Save without XCNV")
    save(gen_cnv(title,ucsc_url,scores,cnv_len,query,cnv_type,exclude_overlaps,morbid_gene_overlaps,decipher_overlaps,clinvar_overlaps,dgv_overlaps,dijon_overlaps,ref_overlaps,None),str_query,organ)
	    
    #XCNV Data
    xcnv_data = xcnv(str_query)
	    
	#Delete temp file
    filename = "cnv/" + str_query
    if organ != "":
        filename += "_" + organ
    conn, bucket = getS3()
    deleteFile(bucket,filename + ".tmp")
    conn.close()
	    
    #Overwrite with XCNV
    print("Save with XCNV")
    save(gen_cnv(title,ucsc_url,scores,cnv_len,query,cnv_type,exclude_overlaps,morbid_gene_overlaps,decipher_overlaps,clinvar_overlaps,dgv_overlaps,dijon_overlaps,ref_overlaps,xcnv_data),str_query,organ)

    #Return data for excel creation
    return {"ref":query.ref,"chromosome":query.chr,"start":query.start,"end":query.end,"type":query.type,"len":query.end-query.start,"final_score":scores["score"],"computed_score":scores["computed_score"],"xcnv":xcnv_data["MVP_score"],"clinvar":scores["clinvar"],"decipher":scores["decipher"],"in_dijon":scores["in_dijon"],"is_piev":scores["is_piev"]}
	

#AWS Lambda Handler
def handler(event, context):
	query = event["headers"]["query"]
	organ = event["headers"]["organ"]
	search(query,organ)
	return 200
	
#Process a BED file, a list of query to return as CSV
def process_bed(username,bed_filepath,ref,organ):
	
	#Get S3
    conn, bucket = getS3()
    data = []
    
    #Get BED file data
    bed_data = readFile(bucket,username + "/_cnv_hub/" + bed_filepath).decode().replace("\r","")
    print(bed_data)
    for l in bed_data.split("\n"):
        if "#" not in l and len(l) > 6:
            d = l.split("\t")
            s = search(ref + ":chr" + d[0] + ":" + d[1] + "-" + d[2] + ":" + d[3] ,organ)
            print(s)
            data.append(s)    
            
    #Write CSV
    csv = "ref;chromosome;start;end;type;len;final_score;computed_score;xcnv;clinvar;decipher;in_dijon;is_piev\n"
    for d in data:
        csv += d["ref"] + ";" + d["chromosome"] + ";" + str(d["start"]) + ";" + str(d["end"]) + ";" + d["type"] + ";" + str(d["len"]) + ";" + str(d["final_score"]) + ";" + str(d["computed_score"]) + ";" + str(d["xcnv"]) + ";" + str(d["clinvar"]) + ";" + str(d["decipher"]) + ";" + str(d["in_dijon"]) + ";" + str(d["is_piev"]) + "\n"
    
    #Save to S3
    setFile(bucket,username + "/_cnv_hub/" + bed_filepath + "." + ref + ".csv",csv)
    conn.close()
    
    #Good return
    return 200
