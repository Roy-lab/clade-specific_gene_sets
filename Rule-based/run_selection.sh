set -u

NEWICK=$1	# newick format tree (txt file)
TARGET=$2	# target node name (e.g. PHYPA, ANC3)
MATRIX=$3	# value matrix file (txt file)
DIRECTION=$4	# "increase" | "decrease"
OUT_CONFIG=$5	# output file name for config file
OUT_MATRIX=$6	# output file name for resultant filtered matrix file


python generate_config.py $NEWICK $TARGET $OUT_CONFIG

python collect_geneset.py $CONFIG $MATRIX $DIRECTION $OUT_MATRIX


