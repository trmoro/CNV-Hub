# UCSC Linker
from cnvannot.common.coordinates import GenomicCoordinates
from cnvannot.annotations.base import base_coordinates_annotation


def ucsc_get_annotation_link(query: GenomicCoordinates) -> dict:
    res = "https://genome.ucsc.edu/cgi-bin/hgTracks?db=" + str(query.ref.replace("xcnv","")) + "&lastVirtModeType=default" \
                                                                            "&lastVirtModeExtraState=&virtModeType" \
                                                                            "=default" \
                                                                            "&virtMode=0&nonVirtPosition=&position=" + \
          str(query.chr) + "%3A" + str(query.start) + "%2D" + str(query.end) + "&hgsid" \
                                                                               "=1321310419_kHoGpB0FxWh71zY7ManBwpAd7NkE "

    res_dict = base_coordinates_annotation(query)
    res_dict["ucsc"] = {"link": res}

    return res_dict
