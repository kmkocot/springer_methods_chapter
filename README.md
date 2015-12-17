README Springer Methods Chapter

Scripts used in Chapter "Phylogenomic approaches using transcriptome data" by J.T. Cannon
and K.M. Kocot. In S. Bourlat, Ed., Marine Genomics – Methods and Protocols. Springer.

These have been tested on Scientific Linux, but aside from potentially necessary changes 
to "rename" commands it should work in any Linux environment (See chapter text for 
discussion of 'rename' command changes in Ubuntu vs. Scientific Linux). 

Absolute paths (where used) may need to be edited.

batch_prep_sequences.sh takes translated amino acid fasta files and ensures that 
headers are appropriate for downstream analysis with HaMStR.

HaMStR_v13_concatenate.sh renames the files output by hamstr into a format appropriate 
for downstream analysis using the phylogenomics_dataset_assembly.sh script. 
Organisms included in the core ortholog set can be added or removed from each OG 
(see end of script) The rename command will need to be changed if you run this script 
on Ubuntu.

phylogenomics_dataset_assembly.sh takes the output of HaMStR and performs several steps 
to remove groups and sequences that are not suitable for phylogenomic analysis. 

The final product of this script is a set of trimmed amino acid alignments representing
putatively orthologous groups suitable for phylogenomic analysis.

A number of programs must also be in the path including Aliscore, Alicut, MAFFT, FastTree
PhyloTreePruner, HaMStR (for nentferner.pl).

A number of variables must be modified for your purposes within the bash script. 
We suggest you examine the entire script carefully and modify it as needed.

Input fasta file headers must be in the following format: 
>orthology_group_ID|species_name_abbreviation|annotation_or_sequence_ID_information

Example: >0001|LGIG|Contig1234

Fasta headers may not include spaces or non-alphanumeric characters except 
for underscores (pipes are OK as field delimiters only).


