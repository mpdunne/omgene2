from typing import Dict, List, Tuple

from omgene.utils.gff import GeneContext


def get_exon_map(gc: GeneContext) -> str:
    """
    Get a map of the exons, as a string representation. Exons are marked as X, non-exonic regions
    are marked with . and capped with \\ and /.

    :param gc: The GeneContext from which to get the exon map.
    :return: A string representation of the exons in the gene region.
    """
    exons = gc.meta_exons
    chrom = gc.meta_chr_seq

    exon_map = ''
    last_exon_end = 0

    for i, exon in enumerate(exons):
        if i == 0 and exon[0] > 0:
            exon_map += '.' * (exon[0] - last_exon_end - 1) + '/'
        else:
            exon_map += '\\' + '.' * (exon[0] - last_exon_end - 2) + '/'

        exon_length = exon[1] - exon[0]
        exon_map += 'X' * exon_length
        last_exon_end = exon[1]

    exon_map += '\\' + (len(chrom) - last_exon_end - 1) * '.'
    return exon_map


def exons_from_exon_map(exon_map: str) -> List[List[int]]:
    """
    Convert an exon map into a list of exon annotations in the format [start, end].

    :param exon_map: A string representation of some exons in a gene region.
    :return: The coordinates of the exons relative to the start of the gene region.
    """
    exons = []
    position = 0
    current_exon = []

    for i, char in enumerate(exon_map):

        if char == 'X' and current_exon == []:
            current_exon.append(position)

        if (char != 'X' or i == len(exon_map) - 1) and current_exon != []:
            current_exon.append(position)
            exons.append(current_exon)
            current_exon = []

        position += 1

    return exons


def exon_map_to_aa_map(exon_map: str) -> Dict[str, List[int]]:
    """
    For each amino acid given by the exon map, list the coordinates of the exon map that specify
    that amino acid.

    :param exon_map: An exon map.
    :return: An AA map in the form {aa_pos: [x1, x2, x3] for aa_pos in ....}
    """
    aa_map = {}
    aa_ix = 0
    cds_len = 0
    for i, x in enumerate(exon_map):
        if x == 'X':
            aa_map[aa_ix] = aa_map.get(aa_ix, []) + [i]
            cds_len += 1
            if cds_len % 3 == 0:
                aa_ix += 1

    return aa_map


def aligned_seq_region_to_exon_map_region(aligned_seq: str, region: List[int], aa_map: Dict[str, List[int]]) \
        -> Tuple[int]:
    """
    Use an AA map to get the corresponding exon map region from a specified aa map region.

    :param aligned_seq: An aligned sequence.
    :param region: The region of the AA map that we're interested in.
    :param aa_map: An AA: [x1, x2, x3] map.
    :return: The coordinates of which parts of the exon map correspond to the specified region in the AA map.
    """
    l, r = region
    unaligned_l = len([x for x in aligned_seq[:l] if x != '-'])
    unaligned_r = len([x for x in aligned_seq[:r] if x != '-'])
    map_l = aa_map[unaligned_l]
    map_r = aa_map[unaligned_r]
    return map_l[0], map_r[-1]


class TransplantError(Exception):
    pass


def exon_map_transplant(patient: str, donor: str, region: List[int]) -> str:
    """
    Transplants the map contents in the specified region from the donor to the patient.
    If the regions are incompatible, throw an error!

    :param patient: The exon map to which we want to transplant the region.
    :param donor: The exon map from which we want to transplant the healthy region.
    :param region: The coordinates of the region to transplant.
    :return: The exon map for the patient with the donor region transplanted.
    """
    if len(patient) != len(donor):
        raise TransplantError('Exon maps must have the same size in order to perform a transplant.')

    s, e = region
    patient_l = patient[:s]
    patient_region = patient[s:e + 1]
    patient_r = patient[e + 1:]

    donor_l = donor[:s]
    donor_region = donor[s:e + 1]
    donor_r = donor[e + 1:]

    patient_exon_content = len([x for x in patient_region if x == 'X'])
    donor_exon_content = len([x for x in donor_region if x == 'X'])

    if (donor_exon_content - patient_exon_content) % 3 != 0:
        raise TransplantError('Cannot transplant gene model part - incompatible frames.')

    bad_transplant_start = (s != 0) and {patient_l[-1], donor_l[-1], donor_region[0]} == {'X', '.'}
    bad_transplant_end = (s != len(patient)) and {patient_r[0], donor_r[0], donor_region[-1]} == {'X', '.'}

    if bad_transplant_start or bad_transplant_end:
        raise TransplantError('Transplant would introduce an unseen splice site or invalid start codon.')

    return patient_l + donor_region + patient_r
