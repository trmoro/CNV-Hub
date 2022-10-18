#Compute scores, return [global,xcnv score,omim score,decipher score,clinvar score]
def compute_scores(cnv_len,dgv_overlaps,clinvar_overlaps,decipher_overlaps,piev_overlaps,xcnv_score,xcnv_on=True):
	
	criteria = { "cmaj1": [], "cmin2" : [], "cmaj2" : [], "cmin4" : [], "cmaj4" : [], "cmaj5" : [] }
	
	#Length > 1 Mb = cmin4
	if cnv_len > 1000000:
		criteria["cmin4"].append("Length is higher than 1 Mb")
		
	#No gene = cmaj2
	if len(clinvar_overlaps) == 0:
		criteria["cmaj2"].append("No Morbid Gene overlaped")
		
	#DGV > 0.8 = cmaj2, DGV > 0.5 = cmin2
	maxDgv = 0
	for d in dgv_overlaps:
		if d["overlaps"] > maxDgv:
			maxDgv = d["overlaps"]
	if maxDgv > 0.8:
		criteria["cmaj2"].append("DGV overlaps more than 80 % the CNV")
	elif maxDgv > 0.5 :
		criteria["cmin2"].append("DGV overlaps more than 50 % the CNV")
		
	#Clinvar score
	clinvar_score = 0
	clinvar_coef = 0
	clinvar_scores = []
	for c in clinvar_overlaps:
		
		#All "if" because 1 overlap can have two patho state (ex : Benign/Likely_benign)
		
		if "Benign" in c["isPatho"]:
			clinvar_score += 0
			clinvar_scores.append(0)
		if "Likely_benign" in c["isPatho"]:
			clinvar_score += 0.25
			clinvar_scores.append(0.25)
		if "Uncertain_significance" in c["isPatho"]:
			clinvar_score += 0.5
			clinvar_scores.append(0.5)
		if "Likely_pathogenic" in c["isPatho"]:
			clinvar_score += 0.75
			clinvar_scores.append(0.75)
		if "Pathogenic" in c["isPatho"]:
			clinvar_score += 1
			clinvar_scores.append(1)
			
	if len(clinvar_scores) > 0:
		clinvar_score /= len(clinvar_scores)
	else:
		clinvar_score = 0.5
	
	#Clinvar criteria
	if clinvar_score > 0.9:
		criteria["cmaj5"].append("Bibliography : Clinvar pathogenicity score is higher than 80 %")
	elif clinvar_score > 0.75:
		criteria["cmaj4"].append("Bibliography : Clinvar pathogenicity score is higher than 50 %")
	elif clinvar_score < 0.25:
		criteria["cmin2"].append("Bibliography : Clinvar benignity score is higher than 50 %")
	elif clinvar_score < 0.1:
		criteria["cmaj1"].append("Bibliography : Clinvar benignity score is higher than 80 %")

	#Decipher score
	decipher_score = 0
	decipher_coef = 0
	decipher_scores = []
	decipher_coefs = []
	for d in decipher_overlaps:
		
		#All "if" because 1 overlap can have two patho state (ex : Benign/Likely_benign)
		
		if "Benign" in d["pathologic"]:
			decipher_score += 0
			decipher_coef += d["overlaps"]
			decipher_scores.append(0)
			decipher_coefs.append(d["overlaps"])
		if "Likely benign" in d["pathologic"]:
			decipher_score += 0.25 * d["overlaps"]
			decipher_coef += d["overlaps"]
			decipher_scores.append(0.25)
			decipher_coefs.append(d["overlaps"])
		if "Uncertain" in d["pathologic"]:
			decipher_score += 0.5 * d["overlaps"]
			decipher_coef += d["overlaps"]
			decipher_scores.append(0.5)
			decipher_coefs.append(d["overlaps"])
		if "Likely pathogenic" in d["pathologic"]:
			decipher_score += 0.75 * d["overlaps"]
			decipher_coef += d["overlaps"]
			decipher_scores.append(0.75)
			decipher_coefs.append(d["overlaps"])
		if "Pathogenic" in d["pathologic"]:
			decipher_score += d["overlaps"]
			decipher_coef += d["overlaps"]
			decipher_scores.append(1)
			decipher_coefs.append(d["overlaps"])
	
	if decipher_coef > 0:		
		decipher_score /= decipher_coef
	else:
		decipher_score = 0.5
	
	#Decipher criteria
	if decipher_score > 0.9:
		criteria["cmaj5"].append("Bibliography : Decipher pathogenicity score is higher than 80 %")
	elif decipher_score > 0.75:
		criteria["cmaj4"].append("Bibliography : Decipher pathogenicity score is higher than 50 %")
	elif decipher_score < 0.25:
		criteria["cmin2"].append("Bibliography : Decipher benignity score is higher than 50 %")
	elif decipher_score < 0.1:
		criteria["cmaj1"].append("Bibliography : Decipher benignity score is higher than 80 %")
		
	#XCNV
	xcnv = -1
	if xcnv_on:
		xcnv = xcnv_score
		
	#CNV Hub Score
	score = 0
	coef = len(criteria["cmaj1"]) + len(criteria["cmin2"]) + len(criteria["cmaj2"]) + len(criteria["cmin4"]) + len(criteria["cmaj4"]) + len(criteria["cmaj5"])
	score += len(criteria["cmaj2"]) * 0.1
	score += len(criteria["cmin2"]) * 0.25
	score += len(criteria["cmin4"]) * 0.75
	score += len(criteria["cmaj4"]) * 0.9
	score += len(criteria["cmaj5"])
	if coef > 0:
		score /= coef
	else:
		score = 0.5
		
	#PIEV
	piev = 0
	if len(piev_overlaps) > 0:
		score = 0.5
		piev = 1
		
	#Return
	return {"score":score,"clinvar":clinvar_score,"decipher":decipher_score,"xcnv":xcnv,"criteria":criteria,"is_piev":piev,"clinvar_scores":clinvar_scores,"decipher_scores":decipher_scores,"decipher_coefs":decipher_coefs}
		
