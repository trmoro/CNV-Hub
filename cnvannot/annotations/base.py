from cnvannot.common.coordinates import GenomicCoordinates


def base_coordinates_annotation(query: GenomicCoordinates) -> dict:
    return {"base": {"ref": query.ref, "chr": query.chr, "start": query.start, "end": query.end}}
