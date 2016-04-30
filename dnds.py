from __future__ import print_function
from codons import codons
from fractions import Fraction

BASES = {'A', 'G', 'T', 'C'}


def split_seq(seq, n=3):
    '''Returns sequence split into chunks of n characters, default is codons'''
    return [seq[i:i+n] for i in range(0, len(seq), n)]


def average_list(l1, l2):
    return [(i1 + i2) / 2 for i1, i2 in zip(l1, l2)]


def dna_to_protein(codon):
    '''Returns single letter amino acid code for given codon'''
    return codons[codon]


def is_synonymous(codon1, codon2):
    '''Returns boolean whether given codons are synonymous'''
    return dna_to_protein(codon1) == dna_to_protein(codon2)


def dnds_codon(codon):
    '''Returns list of synonymous counts'''
    syn_list = []
    for i in range(len(codon)):
        base = codon[i]
        other_bases = BASES - {base}
        syn = 0
        for new_base in other_bases:
            new_codon = codon[:i] + new_base + codon[i + 1:]
            syn += int(is_synonymous(codon, new_codon))
        syn_list.append(Fraction(syn, 3))
    return syn_list


def dnds_codon_pair(codon1, codon2):
    return average_list(dnds_codon(codon1), dnds_codon(codon2))


def syn_sum(seq1, seq2):
    syn = 0
    codon_list1 = split_seq(seq1)
    codon_list2 = split_seq(seq2)
    for i in range(len(codon_list1)):
        codon1 = codon_list1[i]
        codon2 = codon_list2[i]
        dnds_list = dnds_codon_pair(codon1, codon2)
        syn += sum(dnds_list)
    return syn


def dnds(seq1, seq2):
    # Check that both sequences have the same length
    assert len(seq1) == len(seq2)
    # Check that sequences are codons
    assert len(seq1) % 3 == 0
    assert len(seq2) % 3 == 0
    # Strip any whitespace from both strings
    seq1 = seq1.replace(' ', '')
    seq2 = seq2.replace(' ', '')
    syn = syn_sum(seq1, seq2)
    non = len(seq1) - syn


if __name__ == '__main__':
    dnds('ACC GTG GGA TGC ACC GGT GTG CCC',
         'ACA GTG AGA TAT AAA GGA GAG AAC')
