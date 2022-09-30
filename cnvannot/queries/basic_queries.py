from intervaltree import IntervalTree

from cnvannot.common.coordinates import GenomicCoordinates


def overlap_size(sta1, sta2, end1, end2):
    return max(0, min(end1, end2) - max(sta1, sta2) + 1)

def query_overlaps_get(db, query: GenomicCoordinates) -> bool:
    if query.chr in db:
        itree: IntervalTree = db[query.chr]
        if itree.overlaps(query.start, query.end):
            return itree.envelop(query.start,query.end)

#Compute Overlaps with the overlaps rate
def compute_overlaps(db, query: GenomicCoordinates, envelop=True, minRate = 0):
    res = []
    if query.chr in db:
        if db[query.chr].overlaps(query.start, query.end):
            for r in db[query.chr][query.start:query.end]:
                if not 'var_type' in r.data or query.type == r.data['var_type']:
                    if not envelop or (r.begin >= query.start and r.end <= query.end):
                        overlap_rate = overlap_size(r.begin, query.start, r.end, query.end) / (query.end - query.start)
                        if overlap_rate > minRate:
                            r.data['overlaps'] = overlap_rate
                            res.append(r[2])
    return res             
