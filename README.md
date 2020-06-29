Codon Trees
==========
A Python system to build bacterial phylogenies exploiting PATRIC homology groups (PGFams) and genome annotations. It uses RAxML for tree building.

Given a list of genome IDs (in a file, one per line), the code will request from PATRIC all the gene families for thoss genomes (--homologyScope 'global' yields PGFams, 'local yields PLFams).  It will analyse the distribution of homology groups among genomes and select genes that are single-copy (or within specified limits of allowed duplications and absences).

The target number of genes can be specified to limit the run times. Good performance can be achieved with 10 to 20 genes. Resolving close relationships may benefit from requesting 100 or more genes.

It will align the protein sequences for the selected genes using the program mafft (or muscle if specified by the --aligner option).
In the case of multiple gene copies in a single genome, when allowed by duplication tolerances, a single one is selected based on highest similarity to the other genes in the alignment.

It will then align the DNA sequences to the protein alignment, 3 nucleotides per amino acid, yielding nucleotide sequences aligned on a per-codon basis. This assumes that the coding sequences are strictly in frame.

After aligning codons to proteins, protein alignments are end-trimmed to a specified proportion of occupancy (devault 0.5). Any end positions which are more gaps than this are trimmed back. Nucleotide alignments are end-trimmed to match.

The protein and codon alignments are concatenated into a large data matrix for phylogenetic analysis by RAxML, allowing different rates for the amino acids and 1st, 2nd, and 3rd codon positions.

If the --proteinModel option is AUTO (the default), then an initial RAxML analysis is performed on the aligned proteins to find the best protein model.

The 'fasta bootstrapping' option of RAxML can be enabled to put support values on the tree. If --bootstrapReps is set to zero, then 'fast boostrapping' is not performed but branch support values are found using the RELL method (http://www.ncbi.nlm.nih.gov/pubmed/23418397) using the "-f D" option of RAxML.

Requirements
------------
1. Python 2.7 (and modules: requests, urllib, a few others)
2. BioPython
3. RAxML (https://sco.h-its.org/exelixis/web/software/raxml/index.html)
4. Muscle (https://www.drive5.com/muscle/) or mafft (https://mafft.cbrc.jp/alignment/software/).
5. Figtree (http://tree.bio.ed.ac.uk/software/figtree) to automate generation of graphics. (Optional)

Installation
------------
Have phylocode.py and patric_api.py on the $PYTHONPATH environment variable.
Have raxml and mafft or muscle on the $PATH (you can specify an alternate name for raxml, such as raxml-HPC).
Have the figtree executable on the $PATH, or specify the path to the jar file (--pathToFigtreeJar).

Usage
-----
usage: p3x-build-codon-tree.py [-h] [--parametersJson file.json]
                               [--outputBase filebase]
                               [--outputDirectory out_dir]
                               [--genomeIdsFile [file [file ...]]]
                               [--genomeGroupName [name [name ...]]]
                               [--genomeObjectFile file]
                               [--genomePgfamGeneFile file]
                               [--optionalGenomeIdsFile file]
                               [--homolog_scope global/local] [--maxGenes #]
                               [--excessGenesProp prop] [--excessGenesFixed #]
                               [--bootstrapReps #] [--maxGenomesMissing #]
                               [--maxAllowedDups maxDups] [--aligner program]
                               [--endGapTrimThreshold maxPropGaps]
                               [--raxmlExecutable program_name]
                               [--rateModel model] [--proteinModel model]
                               [--analyzeCodons] [--analyzeProteins]
                               [--threads T] [--deferRaxml]
                               [--writePgfamAlignments] [--writePgfamMatrix]
                               [--pathToFigtreeJar path]
                               [--universalRolesFile path] [--focusGenome id]
                               [--debugMode] [--authToken STRING]
                               [--ignoreAuthEnv] [--ignoreAuthRC]

Codon-oriented aligment and tree analysis of PATRIC protein families

optional arguments:
  -h, --help            show this help message and exit
  --parametersJson file.json
                        parameters in json format (command line overrides)
                        (default: None)
  --outputBase filebase
                        base name for output files, def=codontree (default:
                        None)
  --outputDirectory out_dir
                        for output, create if needed (default: None)
  --genomeIdsFile [file [file ...]]
                        file with PATRIC genome IDs, one per line (or first
                        column of TSV) (default: None)
  --genomeGroupName [name [name ...]]
                        name of user's genome group at PATRIC (default: None)
  --genomeObjectFile file
                        genome object (json file) (default: None)
  --genomePgfamGeneFile file
                        read geneIDs per PGFam per genome from this file
                        (default: None)
  --optionalGenomeIdsFile file
                        optional genome ids, one per line (or first column of
                        TSV) (default: None)
  --homolog_scope global/local
                        use PGFams (global) or PLFams (local} (default:
                        global)
  --maxGenes #          number of genes in concatenated alignment (default:
                        50)
  --excessGenesProp prop
                        multiplier of maxGenes to add to filter out low-
                        scoring alignments (default: 0.5)
  --excessGenesFixed #  fixed excess genes to add to filter out low-scoring
                        alignments (default: 20)
  --bootstrapReps #     number of raxml 'fast boostrap' replicates (default:
                        100)
  --maxGenomesMissing #
                        genomes allowed to lack a member of any homolog group
                        (default: 0)
  --maxAllowedDups maxDups
                        duplicated gene occurrences allowed within homolog
                        group (default: 0)
  --aligner program     program to align protein sequences (default: mafft)
  --endGapTrimThreshold maxPropGaps
                        stringency of end-gap trimming, 0-1.0, lower for less
                        trimming (default: 0.5)
  --raxmlExecutable program_name
                        program to call, possibly with path (default: raxml)
  --rateModel model     variable rate category model CAT|GAMMA (default: CAT)
  --proteinModel model  raxml protein substitution model (default: AUTO)
  --analyzeCodons       analyze only codon nucleotides (default: False)
  --analyzeProteins     analyze only amino acids (default: False)
  --threads T, -t T     threads for raxml (default: 2)
  --deferRaxml          does not run raxml (default: False)
  --writePgfamAlignments
                        write fasta alignment per homolog used for tree
                        (default: False)
  --writePgfamMatrix    write table of counts per homolog per genome (default:
                        False)
  --pathToFigtreeJar path
                        not needed if figtree executable on path (default:
                        None)
  --universalRolesFile path
                        path to file with universal roles to select conserved
                        genes (default: None)
  --focusGenome id      to be highlighted in color in Figtree (default: None)
  --debugMode           more output to log file (default: False)
  --authToken STRING    patric authentication token (default: None)
  --ignoreAuthEnv       turn off authorization by environmental variable
                        (default: False)
  --ignoreAuthRC        turn off authorization by file (default: False)


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

The first file, ending in 'genomeGenePgfams.txt', has all the PGFam and PATRIC gene ids (fig|###.peg.###) for each input genome ID.

The second file, ending in singlishCopyPgfams.txt, has the PGFam IDs selected for analysis, those meeting the criteria for maximum missing genomes and maximum duplicated genes, and truncated to the maximum number of genes specified (--maxGenomesMissing, --maxAllowedDups, --maxGenes).

The file ending in "8taxa_10cds_10proteins.pgfamsAndGenesIncludedInAlignment.txt" specifies the genes included in each homology group, separately by protein and by codon alignment. Sometimes a gene fails in the codon alignment process and will be present in the protein alignment but missing in the other.
The "8taxa" in the file name indicates that of the 9 taxa in the input file, one was missing in the data requested from PATRIC. This sometimes happens.

The file ending in "_codonAndProteins.phy" contains the PHLYP-formatted concatenated alignment for all codons and proteins (or just one or the other if specified). 

The file ending in "_codonAndProteins.partitions" tells raxml what the distinct partitions of data are that are to be analyzed separately, namely the 3 codon positions of DNA and the columns of the alignment with the amino acid characters. It also specifies the substitution matrix to be used for the protein data.

DNA, codon1 = 1-11475\3  
DNA, codon2 = 2-11475\3  
DNA, codon3 = 3-11475\3  
LGF, proteins = 11476-15300  

The file ending in ".raxmlCommand.sh" contains the command line for running raxml on the prepared alignment data. It can be run manually if wanted. E.g.:

raxml -s test9_genome_8taxa_10cds_10proteins_codonAndProteins.phy -n test9_genome_8taxa_10cds_10proteins_codonAndProteins -m GTRGAMMA -q test9_genome_8taxa_10cds_10proteins_codonAndProteins.partitions -p 12345 -T 1 -f a -x 12345 -N 100

The files begining with "RAxML_" are the raxml output files. Of these, the one of interest in this example is this:
RAxML_bipartitions.test9_genome_8taxa_10cds_10proteins_codonAndProteins
This is a Newick file including the support values because the --boostrapReps parameter was specified in this case (as 100).
If the program is run without specifying bootstrapping, then this file will not be generated and the Newick file to select is this:
RAxML_bestTree.test9_genome_8taxa_10cds_10proteins_codonAndProteins

Note that the file beginning with "RAxML_bipartitionsBranchLabels." is formatted differently and, for example, will not load into the program FigTree.

Note that the output tree files use as labels the genome IDs. I may add a utility to relabel these to the genome names.
