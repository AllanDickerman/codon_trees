# Phylogenetic Tree Service

## Overview

The bacterial Phylogenetic Tree Service enables construction of custom phylogenetic trees built from user-selected genomes. 
The **Codon Tree** method selects single-copy BV-BRC PGFams and analyzes aligned proteins and coding DNA from single-copy genes using the program RAxML. 
The resulting tree is combined with user-specified metadata on each genome into a phyloxml file which can be viewed in the Archaeopteryx viewer from the user's BV-BRC workspace.
Alternatively, the Newick version of the file can be downloaded and viewed in FigTree or other software.

Given a list of genome IDs, it will analyse the distribution of homology groups (PGFams) among genomes and select gene families that are single-copy (or within user-specified limits of allowed duplications and absences).
In the case of multiple gene copies in a single genome, when allowed by the the duplication tolerance parameter, a single one is selected based on highest similarity to the other genes in the alignment.

The target number of genes can be specified to limit the run times. Good performance can be achieved with 10 to 20 genes. Resolving close relationships may benefit from requesting 100 or more genes.

It will align the protein sequences for single-copy genes using the program mafft and score the alignment for quality. Having evaluated an excess of genes (if available) it will sort the alignments by score and select the requested number of genes favoring the highest alignment scores.

It will then align the DNA sequences to the protein alignments, 3 nucleotides per amino acid, yielding nucleotide sequences aligned on a per-codon basis.

After aligning codons to proteins, protein alignments are end-trimmed to a user-specified occupancy threshold. Any end positions which have more gaps than this are trimmed. Nucleotide alignments are end-trimmed to match.

The protein and codon alignments are concatenated into a large data matrix for phylogenetic analysis by RAxML, allowing different rates for the amino acids and 1st, 2nd, and 3rd codon positions.

A random sample of the aligned proteins is analyzed by RAxML's exhaustive search for the best model of protein substitution, which is then used for the final analysis.

The 'fast bootstrapping' option of RAxML is used to derive support values on the tree.

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


