import itertools
import numpy as np
import re

from Bio.Align import substitution_matrices
from collections import Counter
from scipy.special import binom
from typing import Dict, List, Tuple


class MSAScorer:
    """
    A class for scoring the consistency of MSAs based on BLOSUM matrices.

    """

    def __init__(self, single_gap: float = -1, double_gap: float = 0):
        """
        Constructor for the MSAScorer class.

        :param single_gap: The score for a single gap.
        :param double_gap: The score for a double gap.
        """
        self.aas = 'ACBEDGFIHKMLNQPSRTWVYXZ*-'
        self.blosum_dict = self._make_blosum_dict(single_gap=single_gap, double_gap=double_gap)
        self.blosum_positions = {e: i for i, e in enumerate(self.aas)}
        self.blosum_matrix = self._make_blosum_matrix(self.blosum_dict)

    @staticmethod
    def _make_blosum_dict(single_gap: float = -1, double_gap: float = 0) -> Dict[Tuple[str], float]:
        """
        Returns a fully symmetric version of the blosum matrix, with gap penalties included.

        :param single_gap: (float) The score for a single gap.
        :param double_gap: (float) The score for a double gap.
        """

        blos = substitution_matrices.load('BLOSUM62')
        blos_dict = {}

        for aa1 in blos.alphabet:
            blos_dict[(aa1, '-')] = blos_dict[('-', aa1)] = single_gap

            for aa2 in blos.alphabet:
                blos_dict[(aa2, aa1)] = blos_dict[(aa1, aa2)] = blos[(aa1, aa2)]

        blos_dict[('-', '-')] = double_gap
        return blos_dict

    def _make_blosum_matrix(self, blosdict: Dict[Tuple[str], int]) -> np.ndarray:
        """
        Convert a blosum dict into a blosum matrix.

        :param blosdict: (Dict[Tuple, float]) A dict of substitution distances for pairs of AAs.
        :return: A matrix of pairwise substitution distances.
        """
        blosmat = np.zeros([len(self.aas), len(self.aas)])
        for i, e in enumerate(self.aas):
            for j, f in enumerate(self.aas):
                blosmat[i, j] = blosdict[(e, f)]
        return blosmat

    def col_score(self, column: List[str]) -> float:
        """
        Get the alignment score for an individual column.

        :param column: A column of amino acid values.
        """
        counts = Counter(column)
        column_length = len(column)
        countsmatrix = np.zeros((len(self.aas), len(self.aas)))

        # Else-else
        for aa1, aa2 in itertools.combinations(counts.keys(), 2):
            pos1 = self.blosum_positions[aa1.upper()]
            pos2 = self.blosum_positions[aa2.upper()]
            countsmatrix[pos1, pos2] = counts[aa1] * counts[aa2]

        # Self-self
        for aa in counts:
            pos = self.blosum_positions[aa.upper()]
            countsmatrix[pos, pos] = binom(counts[aa], 2)

        # Don't count things twice.
        scoresmatrix = self.blosum_matrix * countsmatrix
        score = np.sum(scoresmatrix) / binom(column_length, 2)
        return score

    def alignment_score(self, sequences: List[List[str]], omit_empty: bool = False, scaled: bool = False) -> float:
        """
        Return an score for the alignment, calculated by summing column scores
        based on the Blosum matrix.

        :param sequences: (List of SeqRecords) An iterable of SeqRecord items
        :param omit_empty: (bool) Ignore any empty sequences
        :param scaled: (bool) Whether to scale the alignment score by alignment length.
        """
        if omit_empty:
            sequences = [s for s in sequences if not re.match(r'^-+$', s)]

        if not sequences:
            return 0

        seq_len = len(sequences[0])
        if not seq_len:
            return 0

        columns = {i: [s[i] for s in sequences] for i in range(0, seq_len)}
        scores = {i: self.col_score(columns[i]) for i in columns}
        score = sum(scores.values())

        return score if not scaled else score / (1.0 * seq_len)
