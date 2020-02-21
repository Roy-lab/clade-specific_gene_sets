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

./run_selection.sh **data/newic_tree.txt** **Anc3** **increase** **value_matrix.txt** **example_result/Anc3_inc_config.txt** **example_result/Anc3_inc_matrix.txt**


### Requirement: newick_tree.txt file
### Requirement: value_matrix.txt file
