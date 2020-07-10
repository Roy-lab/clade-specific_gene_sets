### run_selection_de_novo.sh
A wrapper shell script for rule-based selection.

### findTransitionGenesets_miss2
Clade-specific gene set selection program

Example usage:
```
./run_selection_de_novo.sh [input_dir] [matrix.txt] [order.txt] [OG.txt] [col_name] [output_dir_name] [threshold] [min_set_size] [max_missing_allow]

./findTransitionGenesets_miss2 [input_dir] [matrix.txt] [order.txt] [OG.txt] [col_name] [output_dir_name] [threshold] [min_set_size] [max_missing_allow]

 - [input_dir]: newick format tree (txt file)
 - [matrix.txt]: value matrix file (txt file)
 - [order.txt]: order of columns of the matrix.txt without Anc columns (txt file)
 - [OG.txt]: Orthogroup (OG) information file (txt file)
 - [col_name]: ID of representative 
 - [ouput_dir_name]: output directory name
 - [threshold]: parameter for hierarchical clustering branch length cutting point, controls number of output clusters
 - [min_set_size]: parameter for minimum size of the output sets (e.g. 5)
 - [max_missing_allow]: parameter for max number of missing allowance
```

*e.g.*
```
./run_selection_de_novo.sh input/ matrix.txt speciesorder.txt OGIDs_all.txt ARATH thr0.3_miss5/ 0.3 5 5

./findTransitionGenesets_miss2 input/ matrix.txt speciesorder.txt OGIDs_all.txt ARATH thr0.3_miss5/ 0.3 5 5
```

### Requirement file 1: matrix.txt (tab delimited)
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

### Requirement file 2: order.txt
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

### Requirement file 3: OG.txt
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
