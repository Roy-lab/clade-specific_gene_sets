### run_findTransitionGenesets.sh
A running shell script for de novo clustering selection.

Example usage:
```
./run_findTransitionGenesets.sh [clustering_threshold] [missing_threshold]
 - clustering_threshold: threshold value for hierarchical clustering, 
   number means branch length of the dendrogram [0, 1], 
   which affects to the assignment of gene sets.
 - missing_threhsold: threshold count number for allowance of missing values in the profile
```

*e.g.*

./run_findTransitionGenesets.sh **0.3** **5**
 - make gene sets by doing hierarchical clustering based on profiles which contain upto 5 missing values (deal with "0"),
 - and assign gene sets by the cluster assignments of branch length at 0.3 of hierarchical clustering.


### Requirements: 
1. findTransitionGenesets_miss
>- C++ program for doing hierarchical clustering by allowing missing cells in the profile.
>- Compile code by executing "make" in the **code/** directory.

2. input/allspecies_clusterassign.txt
>- Cluster assignment matrix. Result file of Arboretum.
>- Note: file name should be "allspecies_clusterassign.txt" since the **findTransitionGenesets_miss** recognize this name.

3. input/speciesorder.txt
>- List of Species names, ordered by species tree (without ancestors)
>- *e.g.*
```
PHYPA
ORYSJ
MAIZE
SOLTU
MEDTR
ARATH
```

4. input/OGID_4_denovo.txt
>- List of orthogroups (OGs), with profile of gene IDs which is ordered by species tree (note that same order with **speciesorder.txt**)
>- If 
