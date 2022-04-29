---------------
### 1. PROGRAM
---------------
### findTransitionGenesets_miss2
Clade-specific gene set selection program

Example usage:
```
./findTransitionGenesets_miss2 [input_dir] [matrix.txt] [order.txt] [OG.txt] [col_name] [output_dir_name] [threshold] [min_set_size] [max_missing_allow]

 - [input_dir]: input directory that contains matrix, order and OG information
 - [matrix.txt]: value matrix file (txt file)
 - [order.txt]: order of columns of the matrix.txt without Anc columns (txt file)
 - [OG.txt]: Orthogroup (OG) information file (txt file)
 - [col_name]: ID of representative 
 - [ouput_dir_name]: output directory name
 - [threshold]: parameter for hierarchical clustering branch length cutting point
 - [min_set_size]: parameter for minimum size of the output sets
 - [max_missing_allow]: parameter for max number of missing allowance
```

*e.g.*
```
./findTransitionGenesets_miss2 input/ matrix.txt speciesorder.txt OGIDs_all.txt ARATH output/ 0.3 5 5
```

### run_selection_de_novo.sh
A wrapper shell script for running findTransitionGenesets_miss2. Assigned the key input arguments.
Look inside the script for the detailed descriptions for the arguments.

Example usage:
```
./run_selection_de_novo.sh [output_dir_name] [threshold] [max_missing_allow]

```

*e.g.*
```
./run_selection_de_novo.sh output/ 0.3 5
```

-------------------
### 2. REQUIREMENTS
-------------------
### Requirement 1: input directory name
>- A directory name that contains matrix, order and OG information (following requirements 2,3,4)

### Requirement 2: matrix.txt file (tab delimited)
>- First row should be legend for data columns
>- format: " **Loci (\t) col1 (\t) col2 (\t) ...** "
>- Cluster ID starts from "0". e.g. if k=3, cluster IDs are 0, 1, 2.
>- **negative value** means **missing value**, 
   which are genes that are existing in orthogroup but don't have measured expression values (-1)
   OR which are not exsiting ortholog in the orthogroup (-2).
>- Refer to input/matrix.txt

 *e.g.*
``` 
Loci	PHYPA	ORYSJ	MAIZE	Anc5	SOLTU	MEDTR	ARATH	Anc4	Anc3	Anc2	Anc1
OG63_1	-2	-2	-2	-2	0	-2	-2	-2	0	0	0
OG63_2	-2	-2	-2	-2	0	-2	-2	-2	0	0	0
OG63_3	-2	-2	-2	-2	-1	-2	-2	-2	0	0	0
OG223_1	2	-2	-2	-2	-2	-2	-2	-2	-2	3	2
OG223_2	-2	3	3	3	-2	-2	-2	3	3	3	2
```

### Requirement 3: order.txt file
>- List of column names without Ancestral columns (Anc#), ordered same as matrix.txt
>- *e.g.*
```
PHYPA
ORYSJ
MAIZE
SOLTU
MEDTR
ARATH
```

### Requirement 4: OG.txt file
>- List of orthogroups (OGs), with profile of gene IDs, ordered same as order.txt
>- format: " **OGID[\t]GeneID,GeneID,GeneID,...** "
>- Duplicated OGs are designated as "_number". *e.g.* OG223_1, OG223_2 
>- Write "NONE" if there's no species gene ID assigned to the OG.
>- *e.g.*
```
OG474_1	PP1S125_69V6,NONE,NONE,NONE,NONE,AT3G13810
OG474_2	PP1S63_51V6,NONE,NONE,PGSC0003DMG400019342,MTR_2g099990,AT5G60470
OG474_3	NONE,Os01g0242200,Zm00001d009030,PGSC0003DMG400024700,MTR_8g017210,AT5G66730
OG474_4	NONE,NONE,Zm00001d039254,PGSC0003DMG400003372,MTR_4g059870,AT5G44160
```

### Requirement 5: representative column name
>- Name of representative column (species) name in the matrix.txt. *e.g.* ARATH
>- OGs that are having expression value in this column will be named as the contents (gene) name of this column (species).


### Requirement 6: output directory name
>- A directory name that will be written out result files.


### Requirement 7: input parameters
>- **threshold**: parameter for hierarchical clustering branch length cutting point, controls number of output cluster sets

>- **min_set_size**: parameter for minimum size number of the output sets (*e.g.* 5), filters out very small cluster sets less than this number

>- **max_missing_allow**: parameter for maximum number of missing allowance, 
   means ***"By including OGs for the clustering, how many missing values are allowed in a OG tree?"***. 
   So, if this value is "0", only OGs that have all genes and values would be used for the clustering and the ohter OGs with at least one missing value will be dropped.
   About missing values in detail, see the description in matrix.txt (requirement 2) above.

-------------
### 3. OUTPUT
-------------
### Output 1: all_genes_clusterassignment_matrix.txt
>- Same as matrix.txt but grouped and segmented by the de novo clustering result.
>- Separated by "DUMMY###" line. *e.g.*
```
DUMMY387        -100    -100    -100    -100    -100    -100    -100    -100    -100    -100    -100
```

### Output 2: all_genesets.txt
>- Cluster results with cluster IDs and contents (gene, OG) IDs.
>- Contents IDs are delimited by "#"
>- *e.g.*
```
Cluster387      AT2G41160#AT1G10900#AT1G08350#OG3072_1#OG15469_1
```

### Output 3: ordered_mean_clusterassign_matrix.txt
>- Cluster mean value matrix, which is a summary for all_genes_clusterassignment_matrix.txt
>- Each clutser sets are represented as a mean value vector of columns (species).
>- Mean value were calculated by taking mean of column-wise (species-wise mean). Negative missing values are not counted.
>- *e.g.*
```
(clusterset387 in all_genes_clusterassignment_matrix.txt)
AT2G41160       2       -1      -1      -1      2       2       4       4       2       2       2
AT1G10900       2       -1      -1      -1      2       4       4       4       2       2       2
AT1G08350       2       -1      -1      -1      2       4       4       4       2       2       2
OG3072_1        2       -1      -1      -1      2       4       -1      4       2       2       2
OG15469_1       -1      -1      -1      -1      2       4       -1      4       2       2       2

(mean vector for clusterset387 in ordered_mean_clusterassign_matrix.txt)
clusterset387   2       -1      -1      -1      2       4       4       4       2       2       2
```
