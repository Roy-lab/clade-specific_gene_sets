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
 1. First step to make config file is to set one specific "objective" clade, which contains one of more ancestral nodes and extant nodes.
 - e.g. clade3 of dicots (clade with Anc3) or clade5 of monocots (clade with Anc5) in the figure.
 2. Then the each contents would be:
 - **target_anc_point**: most top ancestor of the objective clade
 - **object**: all nodes in the objective clade
 - **extant**: most extant nodes in the objective clade
 - **subancestor**: ancestral node within the clade. This is optional, so if "subancestor" doesn't exists, remain it as blank.
 - **ancestor**: ancestral nodes outside/upper of the objective clade
 - **negate**: ancestral and extant nodes outside of the objective clade, which are in the same or underneath level of the objective calde
 3. each node (column) number should be delimited by " "
 4. indicate data file locations
 - **gnlist**: orthogroups list of all 
 - **matrix**: result matrix from Arboretum
   * cluster IDs should correspond to the expression level, i.e. larger number, highly expressed.
   * missing values should be replaced as "0".
