## run_selection.sh
A running shell script for rule-based selection.

Example usage:
```
./run_selection.sh [config_file] [in|de]
 - in: increased case
 - de: decreased case
```

*e.g.*

./run_selection.sh **config_clade1_moss.txt** **in**


## Requirement: config.txt file
This file contains tree information, which is corresponding to the columns of the matrix.

Example of matrix (table at top), species tree with node number corresponding to the column index and config.txt 
![Example](http://pages.discovery.wisc.edu/~jshin/multi-species-proteome/config_making_example.png)
 - each node (column) number should be delimited by " "
 - if "subancestor" doesn't exists, remain it as blank.
