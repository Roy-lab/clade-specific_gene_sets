### run_selection.sh
A running shell script for rule-based selection.

Example usage:
```
./run_selection.sh [newick.txt] [target_node_name] [increase|decrease] [matrix.txt] [ouput_config_name] [output_matrix_name]

 - [newick.txt]: newick format tree (txt file)
 - [target_node_name]: target node name (e.g. PHYPA, ANC3)
 - [increase|decrease]: direction of alteration, input as "increase" or "decrease"
 - [matrixt.txt]: value matrix file (txt file)
 - [ouput_config_name]: output file name for config file
 - [output_matrix_name]: output file name for resultant filtered matrix file
```

*e.g.*
```
./run_selection.sh  data/newick_tree.txt  Anc3  increase  data/value_matrix.txt  example_result/Anc3_inc_config.txt  example_result/Anc3_inc_matrix.txt
```

### Requirement file 1: newick_tree.txt
> Works on generalized types of newick trees.
> * One ancestral node can have only 1 or 2 extant nodes. (not 3 or more children allowed)
> * The file should be prepared as labeled at ancestral node as well.
> * Ancestral nodes should be named as "Anc#".
> * Refer to data/newick_tree.txt
> * Any direction for the ancestral node label is working: AncX(child_A,child_B) or (child_A,child_B)AncX

*e.g.*
``` 
Anc1(PHYPA,Anc2(Anc5(ORYSJ,MAIZE),Anc3(SOLTU,Anc4(MEDTR,ARATH))))
(PHYPA,((ORYSJ,MAIZE)Anc5,(SOLTU,(MEDTR,ARATH)Anc4)Anc3)Anc2)Anc1
```

### Requirement file 2: value_matrix.txt (tab delimited)
> * First row should be legend for data columns: *e.g. OG_NAME (\t) col1 (\t) col2 (\t) ...*
> * Values are meaning expression levels here. The higher value, the higher exression.
> * "0" value means "missing value".
> * Refer to data/value_matrix.txt

 *e.g.*
``` 
OG_NAME	PHYPA	ORYSJ	MAIZE	Anc5	SOLTU	MEDTR	ARATH	Anc4	Anc3	Anc2	Anc1
OG918_2	0	0	0	0	1	6	0	1	1	1	1
OG378_3	1	0	0	0	0	5	0	3	3	3	1
OG455_2	0	2	2	1	4	0	0	0	4	1	1
```

### Result file 1: config.txt
> * Parsed information of newick_tree.txt file with regarding to the target_node that user input
> * see the example_result/Anc3_inc_config.txt

### Result file 2: matrix.txt (tab delimited)
> * Filtered list of OGs from value_matrix.txt based on criteria that user input
> * see the example_result/Anc3_inc_matrix.txt
