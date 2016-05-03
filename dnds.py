from __future__ import print_function, division
from codons import codons
from fractions import Fraction

BASES = {'A', 'G', 'T', 'C'}


def split_seq(seq, n=3):
    '''Returns sequence split into chunks of n characters, default is codons'''
    return [seq[i:i + n] for i in range(0, len(seq), n)]


def average_list(l1, l2):
    """Return the average of two lists"""
    return [(i1 + i2) / 2 for i1, i2 in zip(l1, l2)]


def dna_to_protein(codon):
    '''Returns single letter amino acid code for given codon'''
    return codons[codon]


def translate(seq):
    """Translate a DNA sequence into the 1-letter amino acid sequence"""
    return "".join([dna_to_protein(codon) for codon in split_seq(seq)])


def is_synonymous(codon1, codon2):
    '''Returns boolean whether given codons are synonymous'''
    return dna_to_protein(codon1) == dna_to_protein(codon2)


def dnds_codon(codon):
    '''Returns list of synonymous counts for a single codon.
    Calculations done per the methodology taught in class.
    http://sites.biology.duke.edu/rausher/DNDS.pdf
    '''
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
    """Get the dN/dS for the given codon pair"""
    return average_list(dnds_codon(codon1), dnds_codon(codon2))


def syn_sum(seq1, seq2):
    """Get the sum of synonymous sites from two DNA sequences"""
    syn = 0
    codon_list1 = split_seq(seq1)
    codon_list2 = split_seq(seq2)
    for i in range(len(codon_list1)):
        codon1 = codon_list1[i]
        codon2 = codon_list2[i]
        dnds_list = dnds_codon_pair(codon1, codon2)
        syn += sum(dnds_list)
    return syn


def hamming(s1, s2):
    """Return the hamming distance between 2 DNA sequences"""
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2)) + abs(len(s1) - len(s2))


def codon_subs(codon1, codon2):
    """Returns number of synonymous substitutions in provided codon pair
    Methodology for multiple substitutions from Dr. Swanson, UWashington
    https://faculty.washington.edu/wjs18/dnds.ppt
    """
    diff = hamming(codon1, codon2)
    if diff < 1:
        return 0
    elif diff == 1:
        return int(translate(codon1) == translate(codon2))
    else:
        syn = 0
        for i in range(len(codon1)):
            base1 = codon1[i]
            base2 = codon2[i]
            if base1 != base2:
                new_codon = codon1[:i] + base2 + codon1[i + 1:]
                syn += int(is_synonymous(codon1, new_codon))
                syn += int(is_synonymous(codon2, new_codon))
        return syn / diff


def substitutions(seq1, seq2):
    """Returns number of synonymous and nonsynonymous substitutions"""
    dna_changes = hamming(seq1, seq2)
    codon_list1 = split_seq(seq1)
    codon_list2 = split_seq(seq2)
    syn = 0
    for i in range(len(codon_list1)):
        codon1 = codon_list1[i]
        codon2 = codon_list2[i]
        syn += codon_subs(codon1, codon2)
    return (syn, dna_changes - syn)


def dnds(seq1, seq2):
    """Main function to calculate dN/dS between two DNA sequences"""
    # Strip any whitespace from both strings
    seq1 = seq1.replace(' ', '')
    seq2 = seq2.replace(' ', '')
    # Check that both sequences have the same length
    assert len(seq1) == len(seq2)
    # Check that sequences are codons
    assert len(seq1) % 3 == 0
    assert len(seq2) % 3 == 0
    syn_sites = syn_sum(seq1, seq2)
    non_sites = len(seq1) - syn_sites
    # print(syn_sites, non_sites)
    syn_subs, non_subs = substitutions(seq1, seq2)
    # print('dN: {} / {}\t\tdS: {} / {}'
    #       .format(non_subs, round(non_sites), syn_subs, round(syn_sites)))
    dn = non_subs / non_sites
    ds = syn_subs / syn_sites
    return dn / ds


if __name__ == '__main__':
    print(dnds('ACC GTG GGA TGC ACC GGT GTG CCC',
               'ACA GTG AGA TAT AAA GGA GAG AAC'))
