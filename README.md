# insertions_deletions_RNA
This repository contains the code accompanying the paper "Insertions and deletions in the RNA sequence-structure map" (Nora S. Martin, Sebastian E. Ahnert, 2021).

This code performs the data analysis and plots presented in the paper. All data is calculated using ViennaRNA (2.4.14). The code for analysis and plotting is written in Python 2.7 and the following packages were used: Matplotlib (2.2.3), NetworkX (2.2), numpy (1.16.5), seaborn (0.9.0), pandas (0.24.2) and subprocess (3.5.4).

References:
- Site-scanning is adapted from the algorithm described in Weiß, Marcel; Ahnert, Sebastian E. (2020): Supplementary Information from Using small samples to estimate neutral component size and robustness in the genotype-phenotype map of RNA secondary structure. The Royal Society. Journal contribution. https://doi.org/10.6084/m9.figshare.12200357.v2; here base pair swaps are used as well as point mutations and the method is used for shapes instead of full secondary structures.
- ViennaRNA manual: https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/RNAlib-2.4.14.pdf
- RNAshape definition as described in Figure 2 of Janssen, S., Reeder, J. & Giegerich, R. Shape based indexing for faster search of RNA family databases. BMC Bioinformatics 9, 131 (2008). - hairpin loops were not included in level-1 shapes here.

Data files:
The central file is Psample_seq_struct_dataL30_11_2SL2kbT15_10sample500.csv. This file contains the sequence sample used, the predicted shape for each sequence, its robustness and non-neutral neighbours. The non-neutral neighbours for each sequence are saved in a single field, seperated by 'o' characters.

Nomenclature:
The samples used here are phenotype samples (or "P-samples") in the sense used by Dingle et al. (Dingle, K., Schaper, S., & Louis, A. A. (2015). The structure of the genotype–phenotype map strongly constrains the evolution of non-coding RNA. Interface focus, 5(6), 20150053): each structure (or phenotype) is represented equally since there is an equal number of sequences per neutral set (with the possible exception of small neutral sets, as discussed in the paper). Therefore, the naming conventions of files and folders often starts with a capital "P" for "Psample". Where we work with a (pseudo)random sequence sample, this is denoted as "genotype-sampling" or "G-sampling", also following Dingle et al.



