# usage	: ./run_findTransitionGenesets.sh [hclust_threshold] [missing_threshold]
# e.g.	: ./run_findTransitionGenesets.sh 0.05 0

THR=${1}	# e.g. 0.05
MISSING=${2}	# e.g. 0

ORDER=input/speciesorder.txt
OGID=input/OGID_4_denovo.txt
SRCSPC=ARATH
SIZE=5	# size of cluster (# of contents)

OUTDIR=thr${THR}_miss${MISSING}
mkdir ${OUTDIR}

echo "${OUTDIR}"
cp -r input/ in/
./findTransitionGenesets_miss in/ ${ORDER} ${OGID} ${SRCSPC} ${THR} ${OUTDIR} ${SIZE} ${ORDER} ${MISSING}
mv in/ ${OUTDIR}/

echo " - FIN"
echo ""

