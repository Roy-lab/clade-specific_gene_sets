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
>- C++ program for fulfilling hierarchical clustering by allowing missing cells in the profile.
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
>- List of orthogroups (OGs), with profile of gene IDs which is ordered by species tree 
  (note that same order with **speciesorder.txt**)
>- format: " **OGID[\t]GeneID,GeneID,GeneID,...** "
>- Write "NONE" if there's no species gene ID assigned to the OG.
>- *e.g.*
```
OG474_1	PP1S125_69V6,NONE,NONE,NONE,NONE,AT3G13810
OG474_2	PP1S63_51V6,NONE,NONE,PGSC0003DMG400019342,MTR_2g099990,AT5G60470
OG474_3	NONE,Os01g0242200,Zm00001d009030,PGSC0003DMG400024700,MTR_8g017210,AT5G66730
OG474_4	NONE,NONE,Zm00001d039254,PGSC0003DMG400003372,MTR_4g059870,AT5G44160
```
