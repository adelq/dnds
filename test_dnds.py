from __future__ import division
from fractions import Fraction
from nose.tools import assert_equal, assert_almost_equal
from dnds import dnds, pnps, substitutions, dnds_codon, dnds_codon_pair, syn_sum, translate

# From Canvas practice problem
TEST_SEQ1 = 'ACTCCGAACGGGGCGTTAGAGTTGAAACCCGTTAGA'
TEST_SEQ2 = 'ACGCCGATCGGCGCGATAGGGTTCAAGCTCGTACGA'
# From in-class problem set
TEST_SEQ3 = 'ATGCTTTTGAAATCGATCGTTCGTTCACATCGATGGATC'
TEST_SEQ4 = 'ATGCGTTCGAAGTCGATCGATCGCTCAGATCGATCGATC'
# From http://bioinformatics.cvr.ac.uk/blog/calculating-dnds-for-ngs-datasets/
TEST_SEQ5 = 'ATGAAACCCGGGTTTTAA'
TEST_SEQ6 = 'ATGAAACGCGGCTACTAA'


def test_translate():
    assert_equal(translate(TEST_SEQ1), 'TPNGALELKPVR')
    assert_equal(translate(TEST_SEQ2), 'TPIGAIGFKLVR')
    assert_equal(translate(TEST_SEQ3), 'MLLKSIVRSHRWI')


def test_dnds_codon_easy():
    assert_equal([0, 0, 1], dnds_codon('ACC'))


def test_dnds_codon_harder1():
    assert_equal([0, 0, Fraction(1, 3)], dnds_codon('AAC'))
    assert_equal([0, 0, Fraction(2, 3)], dnds_codon('ATC'))


def test_dnds_codon_harder2():
    assert_equal([Fraction(1, 3), 0, Fraction(1, 3)], dnds_codon('AGA'))
    assert_equal([Fraction(1, 3), 0, 1], dnds_codon('CGA'))


def test_dnds_codon_ngs():
    assert_equal([0, 0, 0], dnds_codon('TGG'))
    assert_equal([Fraction(1, 3), 0, 1], dnds_codon('CGG'))


def test_dnds_codon_ngs_sums():
    assert_equal(sum(dnds_codon('ATG')), 0)
    assert_equal(sum(dnds_codon('AAA')), Fraction(1, 3))
    assert_equal(sum(dnds_codon('CCC')), 1)
    assert_equal(sum(dnds_codon('GGG')), 1)
    assert_equal(sum(dnds_codon('TTT')), Fraction(1, 3))
    assert_equal(sum(dnds_codon('TAA')), Fraction(2, 3))


def test_dnds_codon_ps():
    assert_equal(dnds_codon('CTT'), [0, 0, 1])
    assert_equal(dnds_codon('CGT'), [0, 0, 1])


def test_dnds_codon_pair():
    assert_equal([Fraction(1, 3), 0, Fraction(2, 3)],
                 dnds_codon_pair('AGA', 'CGA'))


def test_dnds_codon_pair_harder1():
    assert_equal([Fraction(1, 6), 0, Fraction(1, 3)],
                 dnds_codon_pair('TTG', 'TTC'))


def test_dnds_codon_pair_harder2():
    assert_equal([Fraction(1, 6), 0, Fraction(1, 2)],
                 dnds_codon_pair('TTA', 'ATA'))


def test_syn_sum_ps():
    assert_almost_equal(syn_sum(TEST_SEQ3, TEST_SEQ4), 10, delta=1)


def test_syn_subs():
    assert_equal(substitutions(TEST_SEQ1, TEST_SEQ2), (5, 5))
    assert_equal(substitutions(TEST_SEQ3, TEST_SEQ4), (2, 5))
    assert_equal(substitutions("ACCGGA", "ACAAGA"), (1, 1))
    assert_equal(substitutions("TGC", "TAT"), (1, 1))
    assert_equal(substitutions("ACC", "AAA"), (0.5, 1.5))
    assert_equal(substitutions("CCC", "AAC"), (0, 2))


def test_pnps():
    assert_almost_equal(pnps(TEST_SEQ1, TEST_SEQ2), 0.269, delta=0.1)
    assert_almost_equal(pnps(TEST_SEQ3, TEST_SEQ4), 0.86, delta=0.1)
    assert_almost_equal(pnps(TEST_SEQ5, TEST_SEQ6), 0.1364 / 0.6001, delta=1e-4)


def test_dnds():
    assert_almost_equal(dnds(TEST_SEQ5, TEST_SEQ6), 0.1247, delta=1e-4)

# DNDS.pdf - wrong based on email conversation with Sean
#
# def test_syn_sum():
#     assert_close(syn_sum(TEST_SEQ1, TEST_SEQ2), 7.5833)
