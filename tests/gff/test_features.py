from omgene.utils.gff import GeneContext
from omgene.utils.features import find_start, find_stop


def test_find_start_left_exists():
    gc = GeneContext([[12, 20, '+']], 'NNNATGNNNNNNAAGTTTTGANNNNNN')
    gc_new = find_start(gc)
    assert gc_new is not None
    assert gc_new.exons == [[3, 20, '+']]


def test_find_start_right_exists():
    gc = GeneContext([[12, 20, '+']], 'NNNATGNNNNNNAAGATGTGANNNNNN')
    gc_new = find_start(gc, direction='right')
    assert gc_new is not None
    assert gc_new.exons == [[15, 20, '+']]


def test_find_start_left_exists_but_with_stop_codon():
    gc = GeneContext([[12, 20, '+']], 'NNNATGTGANNNAAGTTTTGANNNNNN')
    gc_new = find_start(gc)
    assert gc_new is None


def test_find_start_left_does_not_exist():
    gc = GeneContext([[12, 20, '+']], 'NNNAAGNNNNNNAAGTTTTGANNNNNN')
    gc_new = find_start(gc)
    assert gc_new is None


def test_find_start_right_does_not_exist():
    gc = GeneContext([[12, 20, '+']], 'NNNAAGNNNNNNAAGTTTTGANNNNNN')
    gc_new = find_start(gc, direction='left')
    assert gc_new is None


def test_find_stop_exists_last_exon_in_frame():
    gc = GeneContext([[12, 21, '+']], 'NNNAAGNNNNNNAAGTTTTTTNNNTGA')
    gc_new = find_stop(gc)
    assert gc_new is not None
    assert gc_new.exons == [[12, 27, '+']]


def test_find_stop_exists_last_exon_out_of_frame():
    gc = GeneContext([[12, 20, '+']], 'NNNAAGNNNNNNAAGTTTTTTNNNTGA')
    gc_new = find_stop(gc)
    assert gc_new is not None
    assert gc_new.exons == [[12, 27, '+']]


def test_find_stop_does_not_exist():
    gc = GeneContext([[12, 20, '+']], 'NNNAAGNNNNNNAAGTTTTTTNNNNNN')
    gc_new = find_stop(gc)
    assert gc_new is None
