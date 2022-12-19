import re

from typing import List


def get_nonmatching_regions(aln_seq_1: str, aln_seq_2: str, flank_size: int = 5) -> List[List[int]]:
    """
    Find regions with no stretches of flank_size matching characters, flanked by matches of length flank_size

    :param aln_seq_1: The first sequence.
    :param aln_seq_2: The second sequence.
    :param flank_size: The size of the matching flank.
    :return: A list of [start, end] coordinates for any matches.
    """

    if len(aln_seq_1) != len(aln_seq_2):
        raise ValueError('Sequences must be the same length.')

    # Get regions for which the sequences are different

    def represent_matches(s1, s2):
        if s1 == s2 == '-':
            return '-'
        elif s1 == s2:
            return 'M'
        else:
            return '.'

    equality = ''.join(represent_matches(x1, x2) for x1, x2 in zip(aln_seq_1, aln_seq_2))
    equality = 'M' * flank_size + equality + 'M' * flank_size
    equality_reduced = ''.join(
        'm' if all(x == 'M' for x in equality[i:i + flank_size]) else '.' for i in
        range(len(equality) - flank_size))

    matches = re.finditer(fr'm\.+m', equality_reduced)
    nonmatching_regions = []

    for match in matches:
        expanded_coords = [match.start(), match.end() + (flank_size - 1)]
        unflanked_coords = [expanded_coords[0], expanded_coords[1] - flank_size]
        nonmatching_regions.append(unflanked_coords)

    return nonmatching_regions
