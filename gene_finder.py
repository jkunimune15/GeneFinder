# -*- coding: utf-8 -*-
"""
GENE_FINDER ASSIGNMENT

@author: Justin Kunimune

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'

    #Might as well check all four
    >>> get_complement('G')
    'C'
    >>> get_complement('T')
    'A'

    #If it's not a nucleotide, it should throw an error
    SKIP>>> get_complement('triangle')
    'Invalid nucleotide'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'T':
        return 'A'
    else:
        raise Exception('Invalid nucleotide')   # nucleotides must be a,c,g, or t


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'

    #reverse complements should go both ways
    >>> get_reverse_complement("AAAGCGGGCAT")
    'ATGCCCGCTTT'
    >>> get_reverse_complement("TGAACGCGG")
    'CCGCGTTCA'
    """
    rna = ""    # to be the output
    for i in range(1, len(dna)+1):   # cycle through the dna backward
        rna += get_complement(dna[-i])

    return rna


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'

    #the method must not read out-of-phase stop codons
    >>> rest_of_ORF("ATGCTAGGTTGA")
    'ATGCTAGGT'

    #it must transcribe the entirety of incomplete sequences
    >>> rest_of_ORF("ATGCACTA")
    'ATGCACTA'
    """
    output = ""
    for i in range (0,len(dna)/3):
        nextCodon = dna[i*3 : (i+1)*3]
        if (nextCodon == 'TAG' or nextCodon == 'TAA' or nextCodon == 'TGA'):
            return output
        else:
            output += nextCodon
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']

    #must not find out of phase ORFs
    >>> find_all_ORFs_oneframe("ATATGCATCAATTAACGTAGCAATGCTAGA")
    []

    # or nested ORFs
    >>> find_all_ORFs_oneframe("ATGATGTAGTAG")
    ['ATGATG']
    """
    list = []
    i = 0
    while i < len(dna)/3:
        if dna[3*i:3*(i+1)] == 'ATG':
            list.append(rest_of_ORF(dna[i*3:]))
            i = i + len(list[len(list)-1])/3
        i = i + 1

    return list


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("ATGATGTAGTAG") # do not find nested ORFs
    ['ATGATG']
    """
    return find_all_ORFs_oneframe(dna) + find_all_ORFs_oneframe(dna[1:]) + find_all_ORFs_oneframe(dna[2:])


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA") #Included doctest sufficient; it contains both nested and backwards sequences
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    return find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    list = find_all_ORFs_both_strands(dna)
    longest = ""
    for i in range (0, len(list)):
        if len(list[i]) > len(longest):
            longest = list[i]
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
        The simplest example I could think of:
        >>> longest_ORF_noncoding("ATGCCC", 100)
        6
    """
    dnas = [shuffle_string(dna) for i in range(0,num_trials)]
    maxLength = 0
    for dna in dnas:
        length = len(longest_ORF(dna))
        if length > maxLength:
            maxLength = length
    return maxLength


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'

        Provided doctests are sufficient; this program is simple enough
    """
    result = ""
    for i in range(0,len(dna)/3):
        result += aa_table[dna[3*i:3*(i+1)]]
    return result


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
        
        >>> gene_finder("ATGCCC")
        ['MP']
    """
    cutoffLength = longest_ORF_noncoding(dna,1500)/3
    sequences = [coding_strand_to_AA(ORF) for ORF in find_all_ORFs_both_strands(dna)]
    for i in range(len(sequences)-1,-1,-1):     # I had to use a which loop here for it to work.
        if len(sequences[i]) < cutoffLength:    # I don't know why
            sequences.pop(i)
    return sequences


def longest_substring(str1, str2):
    """ Returns the length of the longest common substring between str1 and str2

        str1: a string
        str2: a string
        returns: the length of the longest common stubstring

        The longest s.s. is 'establishment', so the answer is 13
        >>> longest_substring("antidisestablishmentarianism","establishment sanitarianism")
        13

        Longest s.s. is 'TGCATGA', whose length is 7
        >>> longest_substring("ATGCATGCATGA","ATGTTGCATGAC")
        7
    """
    longest_length = 0
    for start1 in range(0,len(str1)):
        for start2 in range(0,len(str2)):   #for every possible pair of starting positions
            length = 0
            while start1+length < len(str1) and start2+length < len(str2) and str1[start1+length] == str2[start2+length]:
                length += 1   # iterate until they stop being the same
            if length > longest_length:  # save the longest lengths
                longest_length = length
    return longest_length # I realize this function could be a lot faster, but I'm too lazy.


def find_nitrogenase(nitrogenase, metagenome):
    """
        Takes a copy of the nitrogenase dna and a list of tuples containing potential matches,
        and returns the id of the best match
        
        >>> find_nitrogenase('NITROGENASE',[('Incfw 10','NITROGENESE'),'Notagene','SPAAAAAAAAACE'])
        'Incfw 10'
    """
    best_length = 0     # the substring length and name of the best match
    best_gene = ""
    for tuple in metagenome:
        similarity = longest_substring(nitrogenase, tuple[1])
        if similarity > best_length:
            best_length = 0
            best_gene = tuple[0]
    return best_gene


if __name__ == "__main__":
    import doctest
    doctest.testmod()

from load import load_seq
dna = load_seq("./data/X73525.fa")
print gene_finder(dna)

#from load import load_nitrogenase_seq
#nitrogenase = load_nitrogenase_seq()
#from load import load_metagenome
#metagenome = load_metagenome()
#print find_nitrogenase(nitrogenase,metagenome)
