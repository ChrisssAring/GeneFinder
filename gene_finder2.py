# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Chris Aring

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
    switcher = {
        "A": "T",
        "G": "C",
        "T": "A",
        "C": "G",
    }
    return switcher.get(nucleotide, "nothing")


def get_reverse_complement(dna):
    splitDNA = list(dna)
    reverseComplement = []
    for DNA in reversed(splitDNA):
        reverseComplement.append(get_complement(DNA))
    return(reverseComplement)


def rest_of_ORF(dna):
    stopCodons = ["TAG", "TGA", "TAA"]
    stopCodonLoc = []
    ORF = []
    codonsInDNA = [dna[i:i+3] for i in range(0, len(dna), 3)]
    for stopCodon in stopCodons:
        stopCodonLoc.extend([i for i, s in enumerate(codonsInDNA)
                             if stopCodon in s])
    stopCodonLoc.sort()
    if not stopCodonLoc:
        return(dna)
    else:
        ORF.append(dna[0:stopCodonLoc[0]])
        return(ORF)


def find_all_ORFs_oneframe(dna):
    index = 0
    allORFs = []
    codonsInDNA = [dna[i:i+3] for i in range(0, len(dna), 3)]
    for codon in codonsInDNA:
        allORFs.extend(rest_of_ORF(dna[index:]))
        index = len(allORFs[codon]) + 3
        print(index)
    return(allORFs)




        if 0 <= i < len(stopCodonLocSorted):
            ORFs.append(dna[(startCodonLocSorted[i]*3+frame):
                            stopCodonLocSorted[i]*3+frame])
        else:
            ORFs.append(dna[(startCodonLocSorted[i]*3+frame):(len(dna))])
        return(ORFs)


def find_all_ORFs(dna):
    frames = [0, 1, 2]
    allORFs = []
    for frame in frames:
        allORFs.append(find_all_ORFs_oneframe(dna, frame))
        return(allORFs)


def find_all_ORFs_both_strands(dna):
    reverseComplement = get_reverse_complement(dna)
    reverseComplementStrandCombined = ''.join(reverseComplement)
    strands = [dna, reverseComplementStrandCombined]
    bothStrands = []
    for strand in strands:
        bothStrands.append(find_all_ORFs(strand))
        return(bothStrands)
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


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
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
