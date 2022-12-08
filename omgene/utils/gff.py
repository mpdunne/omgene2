from typing import Sequence
from Bio.Seq import Seq


class GeneContext:
    """
    A class for containing and manipulating gene coordinates along with the
    sequences of their local genome regions.

    """

    def __init__(self, exons: Sequence[Sequence], chr_seq: str, chr_offset: int = 0):
        """
        Constructor for the GeneContext class.

        :param exons: A list of exons, each in the format (start, end, strand).
        :param chr_seq: The chromosome sequence.
        :param chr_offset: Where this chromosome sequence begins relative to one used in the annotations.
        """

        self._chr_seq = chr_seq
        self._chr_offset = chr_offset

        self.meta_exons = None
        self.meta_chr_seq = None

        self.exons = exons

    def __repr__(self):
        """
        String representation for the GeneContext class.

        :return: String representation.
        """
        all_coords = [x for b, e, s in self.exons for x in (b, e)]
        return f'GeneContext: {len(self.exons)} exons, span {min(all_coords)}-{max(all_coords)}'

    @property
    def seq(self) -> Seq:
        """
        The DNA seq of the annotation.

        :return: The DNA seq of the annotation.
        """
        seq = ''
        for b, e in self.meta_exons:
            seq += self.meta_chr_seq[b: e]
        return Seq(seq)

    @property
    def exons(self):
        """
        The list of exons for the annotation.

        :return: The list of exons for the annotation.
        """
        chr_offset = self._chr_offset
        strand = self.strand
        if strand == '+':
            exons = [[b + chr_offset, e + chr_offset, strand] for b, e in self.meta_exons]
        else:
            chr_len = len(self._chr_seq)
            exons = [[chr_len - e + chr_offset, chr_len - b + chr_offset, strand] for b, e in self.meta_exons][::-1]
        return exons

    @exons.setter
    def exons(self, exons: Sequence[Sequence[str]]):
        """
        Setter for the exons property. Checks that they are consistent and make sense.

        :param exons: The exons to use.
        :return: None.
        """

        # Sort the exons and check that the sequence is monotonic
        exons = sorted(exons, key=lambda x: x[0])
        all_coords = [x for (b, e, s) in exons for x in (b, e)]
        if not all(x <= y for x, y in zip(all_coords, all_coords[1:])):
            raise ValueError('Exons must not overlap.')

        # Check that the exons are in bounds
        too_small = any(x < self._chr_offset for x in all_coords)
        too_big = any(x > self._chr_offset + len(self._chr_seq) for x in all_coords)
        if too_small or too_big:
            raise ValueError('Provided gene model contains out-of-bounds coordinates')

        # Check strand
        strands = [ex[2] for ex in exons]
        if len(set(strands)) != 1:
            raise ValueError('Exons must all be on the same strand.')

        strand = strands[0]
        if strand not in ('-', '+'):
            raise ValueError('Strand must be - or +')
        self._strand = strand

        # Now place the exons into a standard context and save them.
        chr_offset = self._chr_offset
        if strand == '+':
            meta_exons = [[b - chr_offset, e - chr_offset] for b, e, s in exons]
            meta_chr_seq = self._chr_seq
        else:
            chr_len = len(self._chr_seq)
            meta_chr_seq = Seq(self._chr_seq).reverse_complement()
            meta_exons = [[chr_len - e + chr_offset, chr_len - b + chr_offset] for b, e, s in exons][::-1]

        self.meta_exons = meta_exons
        self.meta_chr_seq = meta_chr_seq

    @property
    def strand(self):
        """
        The strand of the annotation.

        :return: The strand of the annotation.
        """
        return self._strand
