---------------
### 1. PROGRAM
---------------
### scFindTransitioning
C++ code for finding transitioning gene sets (hierarchical clustering)

### run_findTransitionGenesets.sh
An examplery wrapper shell script for running scFindTransitioning.<br>
Please see inside the script for the detailed descriptions for the arguments.

Example usage:
```
./run_findTransitionGenesets.sh [input_dir_name] [matrix.txt] [order.txt] [OGID.txt]

```

*e.g.*
```
./run_findTransitionGenesets.sh input/ example_matrix.txt example_order.txt example_OGID.txt
```

-------------------
### 2. REQUIREMENTS
-------------------
### Requirement 1: input directory name
>- A directory name that contains matrix, order and OG information (following requirements 2,3,4)

### Requirement 2: matrix.txt file (tab delimited)
>- Matrix file (text) that contains genes' profiles of patterns (result of Arboretum, usually).
>- Note: Arboretum gives cluster IDs from "0" to "(k-1)". e.g. if run Arboretum with k=3, cluster IDs are 0, 1, 2.
>- 1st row **MUST BE** a legend for data columns
>- 1st row format: " **Loci (\t) c1 (\t) c2 (\t) ...** "
>- Genes' pattern profiles are coming from 2nd row.
>- Format of gene profile row name should have "representative cell cluster ID" and "gene ID", e.g. **c1_ABCC1**
>- profile rows format: " **row_name (\t) pattern_for_c1 (\t) pattern_for_c2 (\t) ...** "
>- Refer to input/matrix.txt

 *e.g.*
``` 
Loci	c1	c3	c4	c2
c1_AAAS	1	1	1	1
c1_ANKRD34B	1	3	3	0
c1_CHRM5	0	3	3	0
c1_CNPY1	0	3	3	0
c1_CRHBP	4	3	4	4
c1_DAB2	1	3	3	0
```

### Requirement 3: order.txt file
>- List of cell cluster names 
>- NOTE: This should be **ordered same as OGID.txt** (Requirement 4)
>- *e.g.*
```
c2
c4
c1
c3
```

### Requirement 4: OG.txt file
>- A legacy input mapping file for indicating the same gene IDs in multiple different contexts (e.g. cell clusters).
>- format: " **OGID[\t]GeneID,GeneID,GeneID,...** "
>- NOTE1: OGID and GeneID are delimited by "\t", GeneIDs in a same row are delimited by ","
>- NOTE2: GeneIDs should be **ordered same as order.txt** (Requirements 3)
>- **OGID**: "OG", line number, followed by "0", e.g. **OG1_0**, **OG42_0**, etc.
>- **GeneID**: "representative cell cluster ID" and "gene ID", e.g. **c1_ABCC1**, **c2_ABCC1**, etc.
>- *e.g.*
```
OG1_0	c2_A1BG,c4_A1BG,c1_A1BG,c3_A1BG
OG2_0	c2_A1BG-AS1,c4_A1BG-AS1,c1_A1BG-AS1,c3_A1BG-AS1
OG3_0	c2_A1CF,c4_A1CF,c1_A1CF,c3_A1CF
OG4_0	c2_A2M,c4_A2M,c1_A2M,c3_A2M
OG5_0	c2_A2M-AS1,c4_A2M-AS1,c1_A2M-AS1,c3_A2M-AS1
```

-------------
### 3. OUTPUT
-------------
>- Output results will be automatically put in a systematically-named directory, such as **"thr0.05_miss0"**.

### Output 1: all_genes_clusterassignment_matrix.txt
>- Same as matrix.txt but grouped and segmented by the de novo clustering result.
>- Separated by "DUMMY###" line. *e.g.*
```
DUMMY387        -100    -100    -100    -100    -100    -100    -100    -100    -100    -100    -100
```

### Output 2: all_genesets.txt
>- Gene sets with transitioning geneset IDs and gene IDs with representative cell cluster (e.g. c1).
>- Gene IDs are delimited by "#"
>- *e.g.*
```
Cluster28	c1_ADARB2#c1_LPIN2#c1_CNTNAP4#c1_PLCH1#c1_PPARGC1A#c1_RIT2#c1_LHFPL2#c1_MDFIC#c1_ZNF385D#c1_AMHR2
```

### Output 3: ordered_mean_clusterassign_matrix.txt
>- A matrix of representative pattern profiles of each transitioning gene set.
>- Patterns are summarized by taking mean values of cluster assignment number of columns (cell clusters)
>- (Note that usually, the cluster assignment number returned from Arboretum is related to the relavant expression level.)
>- *e.g.*
```
Species	c1	c3	c4	c2
clusterset30	1	1	4	1
clusterset32	2	4	4	1
clusterset37	4	4	2	2
clusterset29	2	2	2	4
```
