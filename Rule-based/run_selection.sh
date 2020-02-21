set -u

NEWICK=$1	# newick format tree (txt file)
TARGET=$2	# target node name (e.g. PHYPA, ANC3)
DIRECTION=$3	# "increase" | "decrease"
MATRIX=$4	# value matrix file (txt file)
OUT_CONFIG=$5	# output file name for config file
OUT_MATRIX=$6	# output file name for resultant filtered matrix file


python generate_config.py $NEWICK $TARGET $OUT_CONFIG

python collect_geneset.py $OUT_CONFIG $MATRIX $DIRECTION $OUT_MATRIX
