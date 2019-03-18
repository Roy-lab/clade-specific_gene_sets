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
 - C++ program for doing hierarchical clustering by allowing missing cells in the profile.
 - Compile code by type "make" in the code/ directory.

2. 
