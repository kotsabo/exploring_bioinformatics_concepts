# Bioinformatics 1

The aim of this assignment is to explore some core molecular biology and bioinformatics concepts, and to perform 
manipulations for the investigation of a human disease for which there is a suspected genetic link.

## Question 1
What is the name of the disease you have selected? Explain why it is thought there is a genetic basis for this 
disease. What is the human name for the gene that is thought to be involved? Is this gene known by any other 
names? Whether yes or no, explain how you investigated this. Is there a homologue present in a model system such 
as the mouse, the fruit fly or yeast? If so, give reasons why you think this is a true homologue; if not, 
explain what you did to try to find one.

## Question 2
Now investigate the structure of the gene. Find the gene in a database. Which chromosome is it located on, and 
at which position along the chromosome? What is the structure of the gene, how long is it, how many exons and 
introns does it have, and what percentages are exon and intron? Download the gene transcript and the coding 
sequences for your gene - if multiple are listed, choose the main one, at the top of the list. Where did you 
get this sequence from and what was the unique identifier used so that someone else could be sure they were 
looking at the same sequence? How long is the transcript, and what proportion is coding?

## Question 3
Translate your cDNA sequence into protein/amino acid sequence. How many amino acids does your protein contain? 
Of the 64 possible codons available, how many are used? What is the most common amino acid in the protein? 
How many codons for this amino acid exist and how often is each used?

## Question 4
Now look at the following database:
http://www.kazusa.or.jp/codon/
The codon usage database lists the frequency which each codon is used in a species (different species 
prefer different codons). Sequences which have too many rarer codons result in slowing down transcription 
and inhibition of protein expression - in extreme cases, rare codons are thought to introduce transcription 
errors when the rare tRNA is not available. If you were to try and express your human cDNA sequence in yeast 
(Saccharomyces cerevisiae), which codons in your sequence might cause problems for expression. Note there is 
no hard threshold, but generally codons with 1% usage or less are considered rare.

## Question 5
You are given the following coding sequence fragments. They encode a homologous proteins in different species. 
The sequences are aligned to the correct reading frame:

1. CTGAAGCGGGAGGCTGAGACGCTGCGGGAGCGGGAGGGC
2. CTCAAGCGTGAGGCCGAGACCCTACGGGAGCGGGAAGGC
3. GAAGAGCTGAAGAGAGAGGCTGACAATTTAAAGGACAGA
4. AACGAGGAGCTCAAGCGAGAAGCTGATACGCTGAAGGAC

First, give the likely gene name and the most likely species for these sequences. Sequences 1 and 2 differ 
slightly. How does the resulting protein differ? Could this have functional implications?
Now use the Needleman Wunsch algorithm to compare sequence 1 to the other sequences. Use the scoring: 
match +2, mismatch -3, indel -4. Also perform the comparison with sequence 4 on paper, using the first 
three codons only. Comparing the scores, what would you conclude about the relatedness of the species? 
Is this result consistent with know phylogenetic relatedness?