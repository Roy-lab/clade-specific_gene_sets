set -u

config=$1	# e.g. config_clade1_moss.txt
inORde=$2	# e.g. in / de
script=cladeX_${inORde}creased_selection.pl

./${script} ${config}

