# codon_trees

Codon Trees
==========
A python system to build phylogenies for whole bacterial genomes exploiting PATRIC annotation including PGFams (homology groups).
Given a list of genome IDs (in a file, one per line), the code will request from PATRIC all the PGFams for those genomes.
It will analyse the distribution of PGFams among genomes and select those that are universally single-copy (or within specified limits of duplications and absences).
The maximum number of genes can be specified to limit the run times.
A set of outgroup genomes can be specified which will be included where possible without regard to duplications or absences.
It will align the protein sequences for the selected PGFams using the program muscle.
In the case of multiple gene copies in a single genome, a single one is selected based on highest similarity to the other genes in the alignment.
It will then use the Codon Alignment functionality in BioPython to map the DNA sequences to the protein alignment. 
This yields nucleotide sequences aligned on a per-codon basis.
The codon alignments are concatenated into a large data matrix for phylogenetic analysis by RAxML, allowing different rates for the 
1st, 2nd, and 3rd codon positions.
The protein alignments can also be included in the RAxML data matrix as a distinct partition (and appropriate substitution matrix specified).
The 'fasta bootstrapping' option of RAxML can be enabled to put support values on the tree.
The RAxML output is put into a directory named based on the input file of genome IDs (with "_dir" appended).

Requirements
------------
1. Python 2.7 (and modules: requests, urllib, a few others)
2. BioPython
3. RAxML (https://sco.h-its.org/exelixis/web/software/raxml/index.html)
4. Muscle (https://www.drive5.com/muscle/)

Installation
------------
Have the p3_allan directory with python modules (phylocode.py and patric_api.py)
in a location on the $PYTHONPATH or in the working directory.
Have raxml and muscle on the $PATH (you can specify an alternate name for raxml, such as raxml-HPC)

Usage
-----
usage: buildTreeModular.py [-h] [--outgroupIdsFile file] [--maxGenes #]
                           [--bootstrapReps #] [--maxGenomesMissing #]
                           [--minPropPresent proportion]
                           [--maxAllowedDups maxDups]
                           [--endGapTrimThreshold maxPropGaps]
                           [--raxmlExecutable program_name]
                           [--rateModel model] [--proteinModel model]
                           [--analyzeCodons] [--analyzeProteins]
                           [--analyzeBoth] [--raxmlNumThreads T] [--runRaxml]
                           [--debugMode]
                           genomeIdsFile

positional arguments:
  genomeIdsFile (PATRIC genome ids, one per line, subsequent content on each line ignored)

optional arguments:
  -h, --help            show this help message and exit
  --outgroupIdsFile file
                        ougroup genome ids, one per line (or first column of
                        TSV)
  --maxGenes #          maximum number of genes in concatenated alignment [50]
  --bootstrapReps #     number of raxml 'fast boostrap' replicates [0]
  --maxGenomesMissing #
                        maxumum number of genomes missing a member of the gene family
  --maxGenomesMissing #
                        maximum number of ingroup genomes missing a member of
                        any homolog group
  --maxAllowedDups maxDups
                        maximum duplicated gene occurrences within ingroup
                        genomes for homolog to be included [0]
  --endGapTrimThreshold maxPropGaps
                        stringency of end-gap trimming, lower for less
                        trimming [0.5]
  --raxmlExecutable program_name
                        name of program to call (possibly with path)[raxml]
  --rateModel model     variable rate category model CAT|GTR [CAT]
  --proteinModel model  raxml protein substitution model [WAGF]
  --analyzeCodons       set this flag to analyze codons
  --analyzeProteins     set this flag to analyze proteins
  --analyzeBoth         set this flag to analyze both codons and proteins
  --raxmlNumThreads T   number of threads for raxml [1]
  --runRaxml            set this flag if you want to run raxml (you can run it
                        manually later using the command file provided)
  --debugMode           turns on progress output to stderr
                        
