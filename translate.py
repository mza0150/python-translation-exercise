#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.
    """
    #Crate an empty list to store AA sequence:
    AA_list = []
    # Convert all rna_sequence to upper case:
    rna_sequence=rna_sequence.upper()
    # Convert all rna_sequence into a list:
    rna_list = list(rna_sequence)
    # This conditon will run if rna_sequence is at least 3 bases long, and only once it find start codon ,
    #and stop once it finds stop codon.
    while True:
        if len(rna_list) > 2:
            codon=''.join(rna_list[0:3])
    #Delete first 3 bases since its alread added as codon, thus no longer needed.
            del rna_list[0:3]
        else:
            break
    #Using genetic code dictionary to find AA for each corresponding codon:
        AA=genetic_code[codon]
    #Break loop once it finds stop codon
        if AA=='*':
            break
    #Add add translatable AA to the AA_list:
        AA_list.append(AA)
    return ''.join(AA_list)

def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.
    """
    #pass
    rna_sequence=rna_sequence.upper()
    total_rna_bases=len(rna_sequence)
    last_codon= total_rna_bases-3
    polypeptide_list = []
    for i in range(last_codon+1):
        i_end= i +3
        next_three=rna_sequence[i:i_end]
        if  next_three=='AUG':
            polypeptide=translate_sequence(rna_sequence[i:], genetic_code)
            polypeptide_list.append(polypeptide)
    return polypeptide_list

def get_reverse(sequence):
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, an empty string is returned.
    """
    #pass
    sequence=sequence.upper()
    rna_rev_list=sequence[::-1]
    return rna_rev_list

def get_complement(sequence):
    """Get the complement of `sequence`.

    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, an empty string is returned.
    """
    #pass
    sequence=sequence.upper()
    rna_list=list(sequence)
    comlement_sequence=[]
    complement= {'A' : 'U', 'C' : 'G', 'G': 'C', 'U': 'A'}
    for i in rna_list:
        comlement_sequence.append(complement[i])
    return ''.join(comlement_sequence)

def reverse_and_complement(sequence):
    """Get the reversed and complemented form of `sequence`.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, an empty string is returned.
    """
    #pass
    sequence=sequence.upper()
    rna_list = list(sequence)
    rna_list.reverse()
    rev_c = []
    complement = {'A' : 'U', 'C' : 'G', 'G': 'C', 'U': 'A'}
    for i in rna_list:
        rev_c.append(complement[i])
    return ''.join(rev_c)
def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.
    """
    #pass
    pp_all=[]
    polypeptide_list = get_all_translations(rna_sequence, genetic_code)
    polypeptide_list_long=''.join(polypeptide_list)
    rev_c_seq = reverse_and_complement(rna_sequence)
    polypeptide_list_rev_c=get_all_translations(rev_c_seq, genetic_code)
    polypeptide_list_rev_c_long=''.join(polypeptide_list_rev_c)
    pp_all=polypeptide_list+polypeptide_list_rev_c
    if pp_all==[]:
        return ""
    else:
        return max(pp_all, key=len)

if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq =str("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
