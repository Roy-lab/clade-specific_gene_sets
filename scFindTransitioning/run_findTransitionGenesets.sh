# usage	: ./run_findTransitionGenesets.sh [input_dir_name] [matrix.txt] [order.txt] [OGID.txt]
# e.g.	: ./run_findTransitionGenesets.sh input/ example_matrix.txt example_order.txt example_OGID.txt
set -u

THR=0.05	# Default value is 0.05. Higher threshold value allows the result gene sets having complexed patterns.
MISSING=0	# missing is not allowed in the scFindTransitioning
SIZE=5		# minimum size of cluster (# of genes in a gene set)
SRCCLUST=c1	# Default name of the representative cell cluster

INDIR=$1	# input dir name
MATRIX=$2	# input matrix file, which is the result table of Arboretum
ORDER=$3	# input order file, synced to the OGID file columns
OGID=$4		# input OGID info file

OUTDIR=thr${THR}_miss${MISSING}	# output dir will be named programatically.
mkdir ${OUTDIR}
echo "${OUTDIR}"

./scFindTransitioning ${INDIR} ${MATRIX} ${ORDER} ${OGID} ${SRCCLUST} ${OUTDIR} ${THR} ${SIZE} ${MISSING}
# Meaning of parameters
# ./scFindTransitioning [input_dir] [matrix.txt] [order.txt] [OG.txt] [col_name] [output_dir_name] [threshold] [min_set_size] [max_missing_allow]
#  - [input_dir]: input directory that contains matrix, order and OG information
#  - [matrix.txt]: value matrix file (txt file)
#  - [order.txt]: order of columns of the matrix.txt without Anc columns (txt file)
#  - [OG.txt]: Orthogroup (OG) information file (txt file)
#  - [col_name]: ID of representative 
#  - [ouput_dir_name]: output directory name
#  - [threshold]: parameter for hierarchical clustering branch length cutting point
#  - [min_set_size]: parameter for minimum size of the output sets
#  - [max_missing_allow]: parameter for max number of missing allowance
 
echo " - Clustering finished."
echo ""

