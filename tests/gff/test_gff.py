import pytest

from omgene.utils.gff import GeneContext


def test_gene_context_overlapping_exons():
    exons = [[0, 10, '+'], [8, 16, '+']]
    chr_seq = 'NNNNNNNNNNNNNNNNNNNN'
    with pytest.raises(ValueError) as e:
        _ = GeneContext(exons, chr_seq)
    assert 'Exons must not overlap' in e.value.args[0]


def test_gene_context_exons_out_of_bounds():
    exons = [[0, 10, '+'], [12, 30, '+']]
    chr_seq = 'NNNNNNNNNNNNNNNNNNNN'
    with pytest.raises(ValueError) as e:
        _ = GeneContext(exons, chr_seq)
    assert 'out-of-bounds coordinates' in e.value.args[0]


def test_gene_context_inconsistent_strands():
    exons = [[0, 10, '+'], [12, 16, '-']]
    chr_seq = 'NNNNNNNNNNNNNNNNNNNN'
    with pytest.raises(ValueError) as e:
        _ = GeneContext(exons, chr_seq)
    assert 'Exons must all be on the same strand.' in e.value.args[0]


def test_gene_context_invalid_strands():
    exons = [[0, 10, ';'], [12, 16, ';']]
    chr_seq = 'NNNNNNNNNNNNNNNNNNNN'
    with pytest.raises(ValueError) as e:
        _ = GeneContext(exons, chr_seq)
    assert 'Strand must be - or +' in e.value.args[0]


def test_gene_context_chr_offset_works():
    exons = [[105, 115, '+'], [120, 131, '+']]
    chr_seq = 'NNNNNATGATTTGTCNNNNNATGCTGAACTTNNNNN'
    gc = GeneContext(exons, chr_seq, chr_offset=100)
    assert gc.seq == 'ATGATTTGTCATGCTGAACTT'


@pytest.mark.parametrize("exons,chr_seq,chr_offset", (
        ([[105, 115, '+'], [120, 131, '+']], 'NNNNNATGATTTGTCNNNNNATGCTGAACTTNNNNN', 100),
        ([[105, 116, '-'], [121, 131, '-']], 'NNNNNAAGTTCAGCATNNNNNGACAAATCATNNNNN', 100),
    ))
def test_gene_context_can_recover_exons(exons, chr_seq, chr_offset):
    gc = GeneContext(exons, chr_seq, chr_offset=chr_offset)
    assert gc.meta_exons == [[5, 15], [20, 31]]
    assert gc.exons == exons


@pytest.mark.parametrize("exons,chr_seq,chr_offset", (
        ([[105, 115, '+'], [120, 131, '+']], 'NNNNNATGATTTGTCNNNNNATGCTGAACTTNNNNN', 100),
        ([[105, 116, '-'], [121, 131, '-']], 'NNNNNAAGTTCAGCATNNNNNGACAAATCATNNNNN', 100),
))
def test_gene_context_correct_sequence(exons, chr_seq, chr_offset):
    gc = GeneContext(exons, chr_seq, chr_offset=chr_offset)
    assert gc.meta_exons == [[5, 15], [20, 31]]
    assert str(gc.seq) == 'ATGATTTGTCATGCTGAACTT'
