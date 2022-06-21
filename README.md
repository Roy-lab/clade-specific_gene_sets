# clade-specific_gene_sets (findTransitionGeneset)
The original version with text outputs for Heatmap.awk (http://compbio.mit.edu/pouyak/software/Heatmap.awk).<br><br>
For a sepcificified version of the code, please follow the links:<br>
- de novo method for finding clade specific genesets: https://github.com/Roy-lab/clade-specific_gene_sets/tree/v2/de-novo_clustering
- single cell version (scFindTransitioning): https://github.com/Roy-lab/clade-specific_gene_sets/tree/scFindTransitioning
<br>

## Usage
./findTransitionGenesets_miss [input_dir] [order] [OGID] [SRC_name] [threshold] [output_name] [min_size] [printing_order] [max_missing]

 - [input_dir]: result directory from arboretum/cmint/escarole
 > NOTE for the input cluster assignment matrix file: <br>
 > TAB deimited <br>
 > MUST HAVE HEADER shaped as  " Loci (\t) SPC1 (\t) SPC2 (\t) ... "
 - [order]: order of species (column) names of the matrix.txt WITHOUT Anc columns (a text file)
 - [OGID]: Orthogroup (OG) information file (a text file)
 - [SRC_name]: representative (source) SPC ID name
 - [threshold]: a value for hierarchical clustering branch length cutting point (e.g. 0.05)
 - [ouput_name]: output directory name
 - [min_size]: a value for minimum size of the output sets (e.g. 5)
 - [printing_order]: order of species (column) names for printing out results (a text file)
 - [max_missing]: a values for maximum number of missing allowance (e.g. 0)
<br>

## Specificity of this version (uploaded Jun 21st 2022)
- Does not calling on very big (>100 genes) genesets
- Allow missing values in a gene tree vector
