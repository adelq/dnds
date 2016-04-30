from dnds import dnds_codon, dnds_codon_pair, syn_sum
from fractions import Fraction
from nose.tools import assert_equal, assert_almost_equal

TEST_SEQ1 = 'ACTCCGAACGGGGCGTTAGAGTTGAAACCCGTTAGA'
TEST_SEQ2 = 'ACGCCGATCGGCGCGATAGGGTTCAAGCTCGTACGA'


def test_dnds_codon_easy():
    assert_equal([0, 0, 1], dnds_codon('ACC'))


def test_dnds_codon_harder1():
    assert_equal([0, 0, Fraction(1, 3)], dnds_codon('AAC'))
    assert_equal([0, 0, Fraction(2, 3)], dnds_codon('ATC'))


def test_dnds_codon_harder2():
    assert_equal([Fraction(1, 3), 0, Fraction(1, 3)], dnds_codon('AGA'))
    assert_equal([Fraction(1, 3), 0, 1], dnds_codon('CGA'))


def test_dnds_codon_pair():
    assert_equal([Fraction(1, 3), 0, Fraction(2, 3)],
                 dnds_codon_pair('AGA', 'CGA'))


def test_dnds_codon_pair_harder():
    assert_equal([Fraction(1, 4), 0, Fraction(1, 2)],
                 dnds_codon_pair('TTA', 'ATA'))


def test_syn_sum():
    assert_almost_equal(syn_sum(TEST_SEQ1, TEST_SEQ2), 7.5833)
