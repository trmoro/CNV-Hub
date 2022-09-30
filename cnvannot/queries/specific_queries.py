from cnvannot.common.coordinates import GenomicCoordinates
from cnvannot.queries.basic_queries import overlap_size
from intervaltree import IntervalTree
   
def is_france_incomplete_penetrance(db, query: GenomicCoordinates) -> (bool, str):
    if query.chr in db:
        if db[query.chr].overlaps(query.start, query.end):
            for r in db[query.chr][query.start:query.end]:
                t = r.data['var_type']
                if t == query.type:
                    if overlap_size(r.begin, query.start, r.end, query.end) >= 0.7 * (query.end - query.start):
                        return True, r.data['desc']

    return False, ''


def exc_overlaps_70_percent(db, query: GenomicCoordinates) -> bool:
    if query.chr in db:
        if db[query.chr].overlaps(query.start, query.end):
            for r in db[query.chr][query.start:query.end]:
                if overlap_size(r.begin, query.start, r.end, query.end) >= 0.7 * (query.end - query.start):
                    return True

    return False

#Return the list of all the organs matched in the query
def omim_get_organs(db, query: GenomicCoordinates):
    organList = []
    if query.chr in db:
        itree: IntervalTree = db[query.chr]
        for interv in itree[query.start:query.end]:
            organs = interv.data['organ_list']
            for o in organs.split(";"):
                if o not in organList and o != "":
                    organList.append(o)
    return organList

#Inject Sfari data
def inject_sfari(data, sfari):
    for d in data:
        gene = d["omim_gene_aliases"].split("|")[0]
        if gene in sfari:
            d["sfari"] = sfari[gene]
        else:
            d["sfari"] = {"gene_score" : -1, "syndromic" : -1, "eagle" : -1}

#Inject Gene stat
def inject_gene_stat(data, gene_stat):
    for d in data:
        gene = d["omim_gene_aliases"].split("|")[0]
        if gene in gene_stat:
            d["stats"] = gene_stat[gene]
        else:
            d["stats"] = {"pLI" : -1, "LOEUF" : -1, "sHet" : -1, "pHaplo" : -1, "pTriplo" : -1}
