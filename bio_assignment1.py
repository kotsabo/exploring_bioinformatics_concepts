#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Fri Oct 20 22:38:43 2017

@author: kotsabo
"""

# import Biopython functions
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import matplotlib.pyplot as plt

# load the whole transcipt
transcript = next(SeqIO.parse("sequence.fasta", "fasta"))

# load the coding sequence
coding = next(SeqIO.parse("sequence.txt", "fasta"))

#print(transcript.description)
print("Lenth of transcript: {}.".format(len(transcript.seq)))
print("Coding sequence is {:.4g}% proportion of gene.".format(100*len(coding.seq)/len(transcript.seq)))

#translation of the coding sequence to protein/amino acid sequence
protein = coding.seq.translate()
print(protein)

#number of amino acids
print(len(protein))

protein_codons = list()
different_codons = dict()

for nucleotid in range(0, len(coding.seq), 3):
    codon = str(coding.seq[nucleotid:nucleotid+3])
    protein_codons.append(codon)

    if different_codons.get(codon) != None:
        different_codons[codon] += 1
    else:
        different_codons[codon] = 1

print(len(different_codons))

dict1 = dict(list(different_codons.items())[:len(different_codons)//4])
dict2 = dict(list(different_codons.items())[len(different_codons)//4:len(different_codons)//2])
dict3 = dict(list(different_codons.items())[len(different_codons)//2:3*len(different_codons)//4])
dict4 = dict(list(different_codons.items())[3*len(different_codons)//4:])

def plotBars(d, i):
    plt.clf()
    fig = plt.figure()
    plt.bar(range(len(d)), d.values(), align='center')
    plt.xticks(range(len(d)), d.keys())
    plt.show()
    fig.savefig('fig{}.png'.format(i))
    plt.close()

plotBars(dict1, 1)
plotBars(dict2, 2)
plotBars(dict3, 3)
plotBars(dict4, 4)

amino_acids = dict()

for amino_acid in protein:
    am_ac = str(amino_acid)
    #print(am_ac)
    if amino_acids.get(am_ac) != None:
        amino_acids[am_ac] += 1
    else:
        amino_acids[am_ac] = 1

max_appearance, max_amino_acid = 0, ''

for key, value in amino_acids.items():
    if value > max_appearance:
        max_appearance = value
        max_amino_acid = key
    
print("The amino acid {0} had the maximum appearnaces of {1}.".format(max_amino_acid, max_appearance))

filepath = "results.txt"

def findGene(sequence, filepath):
    file = open("results.txt", "a")
    
    database = 'refseq_rna'
    result_handle = NCBIWWW.qblast("blastn", database, sequence)
        
    # Parse the retuned structure
    blast_records = NCBIXML.parse(result_handle)
    
    # take the first record (we only did one search, so there is only one)
    item = next(blast_records)
        
    E_VALUE_THRESH = 1e-09
    
    for alignment in item.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                output = '****************************\n' \
                        + 'sequence: {0}\n'.format(alignment.title) \
                        + 'length: {0}\n'.format(alignment.length) \
                        + 'e value: {0}\n\n'.format(hsp.expect)
                print(output)
                file.write(output)

    file.close()
                                
findGene('CTGAAGCGGGAGGCTGAGACGCTGCGGGAGCGGGAGGGC', filepath)
findGene('CTCAAGCGTGAGGCCGAGACCCTACGGGAGCGGGAAGGC', filepath)
findGene('GAAGAGCTGAAGAGAGAGGCTGACAATTTAAAGGACAGA', filepath)
findGene('AACGAGGAGCTCAAGCGAGAAGCTGATACGCTGAAGGAC', filepath)





