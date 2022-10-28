import os
import csv
import json

from datetime import datetime

from flask import Flask
from flask_restful import Resource, Api

#Not really used
from cnvannot.annotations.special.france_incomplete_penetrance import france_incomplete_penetrance_load

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
from cnvannot.queries.interpretation import interpretation_get
from cnvannot.queries.specific_queries import exc_overlaps_70_percent, omim_get_organs, inject_sfari, inject_gene_stat
from cnvannot.queries.scoring import compute_scores

#Flask
app = Flask(__name__)
api = Api(app)

france_inc_pen_db = france_incomplete_penetrance_load()

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
        
#XCNV
def do_xcnv(str_query: str):
	
	print("Processing XCNV")
	
	#Process query to .bed
	data = str_query.split(":")
	str_data = data[1].replace("chr","") + "\t" + data[2].replace("-","\t") + "\t" + data[3]
	
	#Write file
	f = open("tmp.bed","w")
	f.write(str_data)
	f.close()
	
	#Execute XCNV and return score
	os.system("./XCNV.R tmp.bed")
	with open("tmp.output.csv") as csvfile:
		csvreader = csv.DictReader(csvfile)
		for row in csvreader:
			return float(row["MVP_score"])
	
	#Remove temp files
	os.remove("tmp.bed")
	os.remove("tmp.sort.bed")
	os.remove("tmp.output.csv")		

#Search
def search(str_query: str, organ: str, xcnv_on = False):
	
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
        
    #XCNV
    xcnv_res = 0
    xcnv_res_interpretation = "not computed"
    if xcnv_on:
        xcnv_res = do_xcnv(str_query)
        xcnv_res_interpretation = xcnv_interpretation_from_score(xcnv_res)
        
	#Compute Scores
    scores = compute_scores(cnv_len,dgv_overlaps,clinvar_overlaps,decipher_overlaps,dijon_overlaps,xcnv_res,xcnv_on)
    #print(scores)
	
    #CNV Synth interpretation
    synth_interpretation = 'Interpretation suggestion(s): ' + interpretation_get(query,
                                                                                 xcnv_res,
                                                                                 exclude_overlaps,
                                                                                 0,
                                                                                 0,
                                                                                 0,
                                                                                 query.type,
                                                                                 france_inc_pen_db,
                                                                                 organ,
                                                                                 organ_list,
                                                                                 0
                                                                                 )

    title = str_query
    if organ != "":
	    title += " (" + organ + ")"
    return {		'title': title,
					'ucsc_url': ucsc_url,
					'score_data': scores,
                    'cnv_len': cnv_len,
                    'cnv_ref': query.ref,
                    'cnv_chr': query.chr,
                    'cnv_start': query.start,
                    'cnv_end': query.end,
                    'cnv_type': cnv_type,
                    'exc_overlaps': exclude_overlaps,
                    'xcnv_res': xcnv_res,
                    'xcnv_res_interpretation': xcnv_res_interpretation,
                    'synth_interpretation': synth_interpretation,
                    'morbid_gene_overlaps': morbid_gene_overlaps,
                    'decipher_overlaps': decipher_overlaps,
                    'clinvar_overlaps': clinvar_overlaps,
                    'dgv_overlaps': dgv_overlaps,
                    'dijon_overlaps': dijon_overlaps,
                    'gene_overlaps': ref_overlaps}

################# API

class cnvSearch(Resource):
	
	def get(self, query):
		
		with open("api.log","a") as logfile:
			logfile.write(datetime.now().strftime("%d/%m/%Y %H:%M:%S") + " : " + query + "\n")
			logfile.close()
		
		res = []
		queries = query.split("&")
		for q in queries:
			organ = ""
			doXcnv = False
			if q.startswith("xcnv"):
				doXcnv = True
				q = q[4:len(q)]
			args = q.split(":")
			if len(args) > 4:
				organ = args[4]
			res.append(search(args[0] + ":" + args[1] + ":" + args[2] + ":" + args[3],organ,doXcnv))
		return res
		
api.add_resource(cnvSearch, "/search/<string:query>")

if __name__ == '__main__':
	app.run(debug=False, host="0.0.0.0")
	
