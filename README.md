# Phylogenetic Tree Service

## Overview

The bacterial Phylogenetic Tree Service enables construction of custom phylogenetic trees built from user-selected genomes. 
The **Codon Tree** method selects single-copy BV-BRC PGFams and analyzes aligned proteins and coding DNA from single-copy genes using the program RAxML. 

Given a list of genome IDs, it will analyse the distribution of homology groups (PGFams) among genomes and select gene families that are single-copy (or within user-specified limits of allowed duplications and absences).
In the case of multiple gene copies in a single genome, when allowed by the the duplication tolerance parameter, a single one is selected based on highest similarity to the other genes in the alignment.

It will align the protein sequences for single-copy genes using the program mafft and score the alignment for quality. Having evaluated an excess of genes (if available) it will sort the alignments by score and select the requested number of genes favoring the highest alignment scores.
The target number of genes can be specified to limit the run times. Good performance can be achieved with 10 to 20 genes. Resolving close relationships may benefit from requesting 100 or more genes.

It will then align the DNA sequences to the protein alignments, 3 nucleotides per amino acid, yielding nucleotide sequences aligned on a per-codon basis.
After aligning codons to proteins, protein alignments are end-trimmed to a user-specified occupancy threshold. Any end positions which have more gaps than this are trimmed. Nucleotide alignments are end-trimmed to match.

The protein and codon alignments are concatenated into a large data matrix for phylogenetic analysis by RAxML, allowing different rates for the amino acids and 1st, 2nd, and 3rd codon positions using the GTR rate model.
An optimal protein substitution model is searched for by RAxML on a sample of the protein alignment columns.
The 'fast bootstrapping' option of RAxML is used to derive support values on the tree.

The resulting tree is combined with user-specified metadata on each genome into a phyloxml file which can be viewed in the Archaeopteryx viewer from the user's BV-BRC workspace.
Alternatively, the Newick version of the file can be downloaded and viewed in FigTree or other software.

A report document is provided describing the genomes, genes and methods used to build the tree. In cases where fewer than the requested number of single-copy genes was available, suggestions are presented as to which genomes could be dropped from the analysis to yield more single-copy genes.

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

This module provides the following application specfication(s):
* [CodonTree](app_specs/CodonTree.md)


## See also

* [Phylogenetic Tree Service Quick Reference](https://www.bv-brc.org/docs/quick_references/services/phylogenetic_tree_building_service.html)
* [Phylogenetic Tree Service](https://www.bv-brc.org/docs/https://bv-brc.org/app/PhylogeneticTree.html)
* [Phylogenetic Tree Service Tutorial](https://www.bv-brc.org/docs//tutorial/phylogenetic_tree/phylogenetic_tree.html)
* [Phylogeny Tab Quick Reference Guide](https://www.bv-brc.org/docs//quick_references/organisms_taxon/phylogeny.html)



## References
1. Davis, J.J., et al., PATtyFams: Protein families for the microbial genomes in the PATRIC database. 2016. 7: p. 118.

2. Han, M.V.; Zmasek, C.M. (2009). "phyloXML: XML for evolutionary biology and comparative genomics". BMC Bioinformatics. 10: 356. doi:10.1186/1471-2105-10-356. PMC 2774328. PMID 19860910

3. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol. 2013 Apr;30(4):772-80. doi: 10.1093/molbev/mst010. 

4. Stamatakis, A., RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 2014. 30(9): p. 1312-1313.

5. Stamatakis, A., P. Hoover, and J. Rougemont, A rapid bootstrap algorithm for the RAxML web servers. Systematic biology, 2008. 57(5): p. 758-771.

6. Zmasek, Christian M.; Eddy, Sean R. (2001). "ATV: display and manipulation of annotated phylogenetic trees". Bioinformatics. 17 (4): 383–384. doi:10.1093/bioinformatics/17.4.383. PMID 11301314

