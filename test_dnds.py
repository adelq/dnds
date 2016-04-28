from dnds import dnds_codon
from nose.tools import assert_equal


def test_dnds_codon():
    assert_equal([0, 0, 1], dnds_codon('ACC'))
