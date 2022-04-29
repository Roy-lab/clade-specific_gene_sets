# usage	: ./run_findTransitionGenesets.sh [hclust_threshold]
# e.g.	: ./run_findTransitionGenesets.sh 0.05

THR=0.05	# Default value is 0.05. Higher threshold value allows the result gene sets having complexed patterns.
MISSING=0	# missing is not allowed in the scFindTransitioning
SIZE=5		# minimum size of cluster (# of genes in a gene set)
SRCCLUST=c1	# Default name of the representative cell cluster

INDIR=input/
MATRIX=example_matrix.txt
ORDER=example_order.txt
OGID=example_OGID.txt

OUTDIR=thr${THR}_miss${MISSING}
mkdir ${OUTDIR}
echo "${OUTDIR}"

cp -r input/ in/
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

