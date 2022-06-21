# clade-specific_gene_sets (findTransitionGeneset)
The original version with text outputs for Heatmap.awk (http://compbio.mit.edu/pouyak/software/Heatmap.awk).<br><br>
For a sepcificified version of the code, please follow the links:<br>
- de novo method for finding clade specific genesets: https://github.com/Roy-lab/clade-specific_gene_sets/tree/v2/de-novo_clustering
- single cell version (scFindTransitioning): https://github.com/Roy-lab/clade-specific_gene_sets/tree/scFindTransitioning
<br>

## Usage
```
./findTransitionGenesets_miss [input_dir] [matrix.txt] [order.txt] [OG.txt] [col_name] [output_dir_name] [threshold] [min_set_size] [max_missing_allow]

 - [input_dir]: result directory from arboretum/escarole
 - [matrix.txt]: value matrix text file (tab deimited, MUST HAVE HEADER shaped as  " Loci (\t) SPC1 (\t) SPC2 (\t) ... ")
 - [order.txt]: order of species (column) names of the matrix.txt WITHOUT Anc columns (a text file)
 - [OG.txt]: Orthogroup (OG) information file (a text file)
 - [col_name]: ID of representative SPC
 - [ouput_dir_name]: output directory name
 - [threshold]: parameter for hierarchical clustering branch length cutting point (e.g. 0.05)
 - [min_set_size]: parameter for minimum size of the output sets (e.g. 5)
 - [max_missing_allow]: parameter for max number of missing allowance (e.g. 0)
```
<br>

## Specificity of this version (uploaded Jun 21st 2022)
- Does not calling on very big (>100 genes) genesets
- Allow missing values in a gene tree vector
