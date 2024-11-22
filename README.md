# DODA_convergence

Sequences, trees, and expression data from Walker-Hale, Guerrero-Rubio & Brockington (_submitted_).

`alignments_and_trees` has the necessary files to reproduce the ancestral sequence reconstruction and convergence results.

`expression` has the necessary data file and script to reproduce normalised fluorescence panels Figure 2.

`structures` has ChimeraX sessions containing AlphaFold2 predicted structures for focal ancestors and representative extant sequences, coloured by branch and lineage substitution type: red, convergent, blue, divergent, dark grey, unique.

- **alignments_and_trees**

  - `Caryophyllales.tre` - the Caryophyllales species topology used for reconciliation.

  - `DODA_full.cds.fa` - FASTA of most inclusive set of DODAa coding sequences.

  - `DODAa_full.cds.aln` - PRANK codon alignment of most inclusive set of DODAa coding sequences.

  - `DODAa_full.cds.aln-cln` - the same alignment after removing columns with < 10% sequence occupancy. This alignment was used for reconciliation.

  - `DODAa_reconciliation.tre` - maximum likelihood reconciliation with GeneRax.

  - `DODAa_reconciliation.pruned.tre` - the same tree after pruning partial sequences.

  - `DODAa_reconciliation.pruned.regrafted.tre` - the same tree after pruning and regrafting sequences from`Stegnosperma`and`Limeum_.

  - `DODAa_pruned.cds.fa` - FASTA of DODAa coding sequences after pruning partial sequences.

  - `DODAa_pruned.cds.aln` - PRANK codon alignment of the same sequences.

  - `DODAa_pruned.pep.aln` - the same alignment translated. This alignment was used for ancestral sequence reconstruction.
  
  - `DODAa_reconciliation.pruned.regrafted.aa_brlen.tre` - the pruned and regrafted reconciliation with amino acid branch lengths used for reconstruction.

  - `DODAa_reconciliation.pruned.regrafted.aa_brlen.node_labels.tre` - the same tree with node labels from FastML.
  
  - `node_correspondence.tsv` - node correspondences for focal sequences in our study. `node_label` gives the FastML node label, `node_id` gives the sequence name, and `node_number` gives the correspondence to the numbers used in the main text. Note that N2 and N5 are separated by a 0-length branch and are identical.

  - `DODAa_pruned.all_reconstructed_nodes.pep.fa` - all reconstructed MAP sequences from FastML. Sequence names correspond to node labels.

  - `DODAa_pruned.focal_reconstructed_nodes.pep.fa` - reconstructed MAP and AltAll sequences for focal nodes. Sequence names correspond to `node_id`.

- **expression**

  - `expression_data.csv` - raw, OD-normalised, and background-corrected fluorescence for yeast strains expressing ancestral and extant sequences.

  - `plotting.R` - script necessary to reproduce panels in Figure 2.

- **structures**

  - `ShDODAa1_7-ShDODAa1_branch_subs.cxs` - predicted structure of ShDODAa1 showing inferred substitutions along the 7-ShDODAa1 branch.

  - `9_8-9_branch_subs.cxs` - predicted structure for node 9 showing inferred substitutions along the 8-9 branch.

  - `10_9-10_branch_subs.cxs` - predicted structure of node 10 showing inferred substitutions along the 9-10 branch.

  - `13_5-13_branch_subs.cxs` - predicted structure of node 13 showing inferred substitutions along the 5-13 branch.

  - `14_13-14_branch_subs.cxs` - predicted structure of node 14 showing inferred substitutions along the 13-14 branch.

  - `15_14-15_branch_subs.cxs` - predicted structure of node 15 showing inferred substitutions along the 14-15 branch.

  - `16_13-16_branch_subs.cxs` - predicted structure of node 16 showing inferred substitutions along the 13-16 branch.

  - `ShDODAa1_7-ShDODAa1_lineage_subs.cxs` - predicted structure of ShDODAa1 showing inferred substitutions along the 7-ShDODAa1 branch (compared to other lineages).

  - `BvDODAa1_8-10_lineage_subs.cxs` - predicted structure of BvDODAa1 showing inferred substitutions along the 8-10 lineage.

  - `McDODAa1_5-15_lineage_subs.cxs` - predicted structure of McDODAa1 showing inferred substitutions along the 5-15 lineage.

  - `CgDODAa1_5-16_lineage_subs.cxs` - predicted structure of CgDODAa1 showing inferred substitutions along the 5-16 lineage.
