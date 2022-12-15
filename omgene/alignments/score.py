import itertools
import numpy as np
import re

from Bio.Seq import Seq
from Bio.Align import substitution_matrices
from collections import Counter
from scipy.special import binom
from typing import Dict, List, Tuple, Union


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

    def _make_blosum_matrix(self, blosdict: Dict[Tuple[str], float]) -> np.ndarray:
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

    def column_scores(self, sequences: List[str]) -> List[float]:
        """
        Return a list of column scores, one for each column.

        :param sequences: (List of SeqRecords) An iterable of SeqRecord items
        """
        columns = [*zip(*sequences)]
        scores = [self.col_score(c) for c in columns]
        return scores

    def alignment_score(self,
                        sequences: List[Union[str, Seq]],
                        omit_empty: bool = False,
                        scaled: bool = False) -> float:
        """
        Return a score for the alignment, calculated by summing column scores
        based on the Blosum matrix.

        :param sequences: (List of SeqRecords) An iterable of SeqRecord items
        :param omit_empty: (bool) Ignore any empty sequences
        :param scaled: (bool) Whether to scale the alignment score by alignment length.
        """
        if len(sequences) == 0:
            return 0.0

        if len(set([len(s) for s in sequences])) != 1:
            raise ValueError('All aligned sequences must have the same length.')

        if omit_empty:
            sequences = [s for s in sequences if not re.match(r'^-+$', str(s))]

        scores = self.column_scores(sequences)
        score = sum(scores)
        seq_len = len(sequences[0])

        return score if not scaled else score / (1.0 * seq_len)
