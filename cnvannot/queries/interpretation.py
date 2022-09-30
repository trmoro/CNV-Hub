from cnvannot.annotations.xcnv import xcnv_interpretation_from_score
from cnvannot.common.coordinates import GenomicCoordinates
from cnvannot.queries.specific_queries import is_france_incomplete_penetrance


def interpretation_get(query: GenomicCoordinates, xcnv_score: float, exc_over: bool, gene_count: int,
                       omim_gene_count: int, dgv_gold_count: int, cnv_type: str, france_db, 
                       given_organ, organ_list, organ_match_count) -> str:
    res = ''

    warning_msg = "\n***The machine interpretation doesn't replace Human interpretation***\n"

    fr_inc_pen = is_france_incomplete_penetrance(france_db, query)
    if fr_inc_pen[0]:
        return fr_inc_pen[1]

    if dgv_gold_count > 0:
        res += 'DGV Databse : Since this variant overlaps one found in the DGV database, the variant could be considered benign\n'

    if organ_match_count > 0:
        res += "Variant overlaps OMIM genes: LIKELY PATHOGENIC"
        if given_organ != "":
            res += " !\n"
        else:
            res += " on "
            for o in organ_list:
                res += o + ", "
            res = res[0:len(res)-2] + " !"
    res += '' if exc_over is False else '\nThis variant overlaps an "Exclusion Region"' '\n' \
                                        'You may consider excluding it'
                                        
    if xcnv_score >= 0.5 and omim_gene_count == 0:
        res += 'Even if X-CNV score is high, meaning that this variant could be pathogenic,\n' \
               'no overlaps with known morbidity-associated genes have been reported,' '\n' \
               'hence, this variant could be benign'
               
    if res == '':
        res = xcnv_interpretation_from_score(xcnv_score).upper()

    #res += warning_msg

    return res
