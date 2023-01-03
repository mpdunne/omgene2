import numpy as np

from copy import deepcopy
from typing import Sequence, Union

from omgene.utils.gff import GeneContext


def find_stop(gene_context: GeneContext, exon_ix: int = -1,
              stop_codons: Sequence[str] = ('TAG', 'TGA', 'TAA')) -> Union[GeneContext, None]:
    """
    Extend the specified exon of the provided gene_context to find a stop codon. Return None if we can't find one.
    TODO: This needs tests. Lots of tests.

    :param gene_context: The gene context.
    :param exon_ix: The index of the exon to extend.
    :param stop_codons: A list of valid stop codons.
    :return: A modified GeneContext with an exon extended into a stop codon.
    """
    gc = deepcopy(gene_context)

    exons = gc.meta_exons
    if abs(exon_ix) > len(exons):
        raise ValueError('Invalid value for exon_ix')

    if exon_ix < 0:
        exon_ix = len(exons) + exon_ix

    exons = gc.meta_exons[:exon_ix + 1]
    last_exon = exons[exon_ix]
    last_exon_end_frame = np.sum([exon[1] - exon[0] for exon in exons]) % 3
    last_exon[1] = int(last_exon[1] - last_exon_end_frame)
    if last_exon[1] < last_exon[0]:
        last_exon[1] += 3

    while True:
        final_codon = gc.meta_chr_seq[last_exon[1] - 3: last_exon[1]]
        if final_codon.upper() in stop_codons:
            exons[-1] = last_exon
            gc.meta_exons = exons
            return gc
        else:
            last_exon[-1] += 3
            if last_exon[-1] >= len(gc.meta_chr_seq):
                return


def find_start(gene_context: GeneContext, exon_ix: int = 0, start_codons: Sequence[str] = ('ATG',),
               stop_codons=('TAG', 'TGA', 'TAA'), direction='left') -> Union[GeneContext, None]:
    """
    Extend the specified exon of the provided gene_context to find a start codon. Return None if we can't find one.
    TODO: This needs tests. Lots of tests.

    :param gene_context: The gene context.
    :param exon_ix: The index of the exon to extend.
    :param start_codons: A list of valid start codons.
    :param stop_codons: Which stop codons to look out for.
    :param direction: Which direction to look in (left or right).
    :return: A modified GeneContext with an exon extended into a start codon.
    """
    gc = deepcopy(gene_context)

    exons = gc.meta_exons
    if abs(exon_ix) > len(exons):
        raise ValueError('Invalid value for exon_ix')

    if exon_ix < 0:
        exon_ix = len(exons) + exon_ix

    first_exon_start_frame = np.sum([exon[1] - exon[0] for exon in exons[:exon_ix]]) % 3
    exons = gc.meta_exons[exon_ix:]
    first_exon = exons[exon_ix]
    first_exon[0] = int(first_exon[0] - first_exon_start_frame)

    if first_exon[1] < first_exon[0] + 3:
        if direction == 'left':
            first_exon[0] -= 3
        else:
            return

    while True:

        first_codon = gc.meta_chr_seq[first_exon[0]: first_exon[0] + 3]
        if first_codon.upper() in start_codons:
            exons[0] = first_exon
            gc.meta_exons = exons
            return gc

        else:
            if first_codon.upper() in stop_codons:
                return

            if direction == 'left':
                first_exon[0] -= 3
            else:
                first_exon[0] += 3

            if first_exon[0] < 0 or first_exon[0] > first_exon[1]:
                return
