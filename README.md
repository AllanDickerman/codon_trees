Codon Trees
==========
A Python system to build phylogenies for whole bacterial genomes exploiting PATRIC annotation including PGFams (homology groups).

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
```
usage: p3x-build-codon-tree.py [-h] [--genomeObjectFile file]
                               [--outgroupIdsFile file] [--maxGenes #]
                               [--bootstrapReps #] [--maxGenomesMissing #]
                               [--maxAllowedDups maxDups]
                               [--endGapTrimThreshold maxPropGaps]
                               [--raxmlExecutable program_name]
                               [--rateModel rateModel]
                               [--proteinModel substModel] [--analyzeCodons]
                               [--analyzeProteins] [--numThreads T]
                               [--deferRaxml] [--outputDirectory out_dir]
                               [--pathToFigtree jar_file]
                               [--focusGenome genome_id] [--debugMode]
                               genomeIdsFile

positional arguments:
  genomeIdsFile         file with PATRIC genome IDs, one per line, optional
                        content after tab delimiter ignored

optional arguments:
  -h, --help            show this help message and exit
  --genomeObjectFile file
                        genome object (json file) to be added to ingroup
                        (default: None)
  --outgroupIdsFile file
                        ougroup genome ids, one per line (or first column of
                        TSV) (default: None)
  --maxGenes #          number of genes in concatenated alignment (default:
                        50)
  --bootstrapReps #     number of raxml 'fast boostrap' replicates (default:
                        0)
  --maxGenomesMissing #
                        ingroup genomes allowed to lack a member of any
                        homolog group (default: 0)
  --maxAllowedDups maxDups
                        duplicated gene occurrences allowed within homolog
                        group (default: 0)
  --endGapTrimThreshold maxPropGaps
                        stringency of end-gap trimming, lower for less
                        trimming (default: 0.5)
  --raxmlExecutable program_name
                        program to call, possibly with path (default: raxml)
  --rateModel rateModel
                        variable rate category model CAT|GAMMA (default: CAT)
  --proteinModel substModel
                        raxml protein substitution model (default: WAGF)
  --analyzeCodons       analyze only codons (ignore amino acids) (default:
                        False)
  --analyzeProteins     analyze only amino acids (default: False)
  --numThreads T        number of threads for raxml (default: 2)
  --deferRaxml          set this flag if you do not want raxml to be run
                        automatically (you can run it manually later using the
                        command file provided) (default: False)
  --outputDirectory out_dir
                        directory for output, create if it does not exist
                        (default: None)
  --pathToFigtree jar_file
                        specify this to generate PDF graphic: java -jar
                        pathToFigtree -graphic PDF CodonTree.nex CodonTree.pdf
                        (default: None)
  --focusGenome genome_id
                        genome to be highlighted in color in Figtree (default:
                        None)
  --debugMode           turns on more progress output to log file (default:
                        False)
```

Input
-----
The only required argument to buildTreeModular.py is the name of a file containing PATRIC genome ids.
There can be a header or not.
These should be one per line, and can be tab-separated from subsequent fields which will be ignored.
E.g.:  
genome.genome_id  
1075089.3  
1171377.3  
1222034.3  
1232659.5  
1249526.3  

Here is the command line used to generate the example output described below:  
python ~/python/buildTreeModular.py --maxGenes 10 --bootstrapReps 100 --maxGenomesMissing 1 --maxAllowedDups 1 --rateModel GAMMA --proteinModel LGF --runRaxml test9_genome.ids 

Output
------
A directory will be created, named by appending "_dir" to the end of the name of the file with the genome IDs (after stripping off a terminal ".ids" if it exists). In the example here, the output directory is named "test9_genome_dir".
All output will go into that directory.  
Here are the files generated by small example based on an input file name of "test9_genomes.id":

```
929K	test9_genome.genomeGenePgfams.txt  
130B	test9_genome.singlishCopyPgfams.txt  
4.0K	test9_genome_8taxa_10cds_10proteins.pgfamsAndGenesIncludedInAlignment.txt  
135K	test9_genome_8taxa_10cds_10proteins_codonAndProteins.phy  
100B	test9_genome_8taxa_10cds_10proteins_codonAndProteins.partitions  
236B	test9_genome_8taxa_10cds_10proteins_codonAndProteins.raxmlCommand.sh  
10K	RAxML_bootstrap.test9_genome_8taxa_10cds_10proteins_codonAndProteins  
20K	RAxML_info.test9_genome_8taxa_10cds_10proteins_codonAndProteins  
480B	RAxML_bipartitionsBranchLabels.test9_genome_8taxa_10cds_10proteins_codonAndProteins  
468B	RAxML_bipartitions.test9_genome_8taxa_10cds_10proteins_codonAndProteins  
454B	RAxML_bestTree.test9_genome_8taxa_10cds_10proteins_codonAndProteins  
```

The first file, ending in 'genomeGenePgfams.txt', has all the PGFam and PATRIC gene ids (fig|###.peg.###) for each input genome ID.

The second file, ending in singlishCopyPgfams.txt, has the PGFam IDs selected for analysis, those meeting the criteria for maximum missing genomes and maximum duplicated genes, and truncated to the maximum number of genes specified (--maxGenomesMissing, --maxAllowedDups, --maxGenes).

The file ending in "8taxa_10cds_10proteins.pgfamsAndGenesIncludedInAlignment.txt" specifies the genes included in each homology group, separately by protein and by codon alignment. Sometimes a gene fails in the codon alignment process and will be present in the protein alignment but missing in the other.
The "8taxa" in the file name indicates that of the 9 taxa in the input file, one was missing in the data requested from PATRIC. This sometimes happens.

The file ending in "_codonAndProteins.phy" contains the PHLYP-formatted concatenated alignment for all codons and proteins (or just one or the other if specified). 

The file ending in "_codonAndProteins.partitions" tells raxml what the distinct partitions of data are that are to be analyzed separately, namely the 3 codon positions of DNA and the columns of the alignment with the amino acid characters. It also specifies the substitution matrix to be used for the protein data.

```
DNA, codon1 = 1-11475\3  
DNA, codon2 = 2-11475\3  
DNA, codon3 = 3-11475\3  
LGF, proteins = 11476-15300  
```

The file ending in ".raxmlCommand.sh" contains the command line for running raxml on the prepared alignment data. It can be run manually if wanted. E.g.:

```
raxml -s test9_genome_8taxa_10cds_10proteins_codonAndProteins.phy -n test9_genome_8taxa_10cds_10proteins_codonAndProteins -m GTRGAMMA -q test9_genome_8taxa_10cds_10proteins_codonAndProteins.partitions -p 12345 -T 1 -f a -x 12345 -N 100
```

The files begining with "RAxML_" are the raxml output files. Of these, the one of interest in this example is this:
RAxML_bipartitions.test9_genome_8taxa_10cds_10proteins_codonAndProteins
This is a Newick file including the support values because the --boostrapReps parameter was specified in this case (as 100).
If the program is run without specifying bootstrapping, then this file will not be generated and the Newick file to select is this:
RAxML_bestTree.test9_genome_8taxa_10cds_10proteins_codonAndProteins

Note that the file beginning with "RAxML_bipartitionsBranchLabels." is formatted differently and, for example, will not load into the program FigTree.

Note that the output tree files use as labels the genome IDs. I may add a utility to relabel these to the genome names.
