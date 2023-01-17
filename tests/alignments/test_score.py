import numpy as np
import pytest
import itertools

from Bio.Align import substitution_matrices
from omgene.alignments.score import MSAScorer


def test_alignment_score_empty_alignment():
    m = MSAScorer()
    seqs = []
    score = m.alignment_score(seqs)
    assert score == 0


def test_alignment_score_only_one_sequence():
    m = MSAScorer()
    seqs = ['AAAAACCCCC']
    score = m.alignment_score(seqs)
    assert score == 0


def test_alignment_score_omit_empty_more_than_two_left():
    m = MSAScorer()
    seqs = ['AAAAACCCCC', 'AAAAACCCCC', '----------']
    score = m.alignment_score(seqs, omit_empty=True)
    assert score != 0


def test_alignment_score_omit_empty_less_than_two_left():
    m = MSAScorer()
    seqs = ['AAAAACCCCC', '----------', '----------']
    score = m.alignment_score(seqs, omit_empty=True)
    assert score == 0


def test_alignment_score_dont_omit_empty():
    m = MSAScorer()
    seqs = ['AAAAACCCCC', '----------', '----------']
    score = m.alignment_score(seqs, omit_empty=False)
    assert np.isclose(score, -6.666666666)


def test_alignment_score_uneven_lengths():
    m = MSAScorer()
    seqs = ['AAAAACCCCC', 'AAAAA']
    with pytest.raises(ValueError):
        _ = m.alignment_score(seqs, omit_empty=False)


def test_alignment_score_correct_value_unscaled():
    m = MSAScorer()
    seqs = ['AAAAACCCCC', 'ACACACACAC']
    score = m.alignment_score(seqs, scaled=False)
    assert score == 39


def test_alignment_score_correct_value_scaled():
    m = MSAScorer()
    seqs = ['AAAAACCCCC', 'ACACACACAC']
    score = m.alignment_score(seqs, scaled=True)
    assert score == 3.9


@pytest.mark.parametrize("column", (
        ('A', 'A', 'A', 'A', 'A'),
        ('A', 'A', 'A', 'C', 'C'),
        ('A', 'C', 'D', 'E', 'F'),
        ('A', 'C', 'D', '-', '-'),
        ('A', 'A', 'A', '-', '-'),
        ('A', 'A', 'A', 'C', 'C', 'C', 'D', 'D', 'E', 'F', '-', '-'),
        ('A', 'A', 'A', 'C', 'C', 'C', 'D', 'D', 'E', 'F', '-', '-'),
    ))
def test_col_score_is_mean_of_pairwise_scores(column):
    blos = substitution_matrices.load('BLOSUM62')
    single_gap_penalty = -1
    double_gap_penalty = 0

    m = MSAScorer()
    col_score = m.col_score(column)
    pairwise_scores = [double_gap_penalty if aa1 == aa2 == '-' else
                       single_gap_penalty if '-' in (aa1, aa2) else
                       blos[aa1, aa2] for aa1, aa2 in itertools.combinations(column, 2)]
    manual_col_score = np.mean(pairwise_scores)
    assert col_score == np.mean(manual_col_score)

