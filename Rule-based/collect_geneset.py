"""
Filtering valued matrix with rules based on the config file.

USAGE:
python collect_geneset.py [config_file] [matrix] [increase|decrease] [output_file_name]

Original: junha.shin@wisc.edu (Feb 2020)
"""
import sys, re

### INPUT ###
config_file = sys.argv[1]	# config file
matrix_file = sys.argv[2]	# matrix
tendency_direction = sys.argv[3]# "increase" | "decrease"
if not tendency_direction == "increase" and not tendency_direction == "decrease":
	sys.exit("ERROR incorrect tendency direction: %s\n - Use only \"increase\" OR \"decrease\"\n" % tendency_direction)
outfile = sys.argv[4]

distinction_level = 2		# level integer of increased or decreased (default=2)
allowed_noise_level = 1		# level integer of noise allowance (default=1)


### LOAD CONFIG DATA ###
config = dict()
f = open (config_file, 'r')
for line in f:
	line = line.rstrip()
	temp = line.split("=")
	config[temp[0]] = temp[1]
f.close()
#print (config)

# Assign vars
target_point_col = config["target_point"]
objects = config["object"].split(" ")
ancestors = config["ancestor"].split(" ")
subancestors = []
if len(config["subancestor"].split(" ")) > 0:
	subancestors = config["subancestor"].split(" ")
negates = config["negate"].split(" ")


### LOAD MATRIX INDEX ###
index = dict()
f = open (matrix_file, 'r')
line1 = f.readline()
f.close()
line1 = line1.rstrip()
#temp1 = np.array(line1.split("\t"))
temp1 = line1.split("\t")
colidx = 0
for col in temp1:
	#index[col] = np.where(temp1 == col)[0][0]
	index[col] = colidx
	colidx += 1
#print (index)



### FILTERING ###
result = []
f = open (matrix_file, 'r')
for line in f:
	line = line.rstrip()
	if re.search("^OG_NAME", line):
		result.append(line)
		continue

	temp = line.split("\t")

	## RULE1: target must have value (!=0), any of negates must have value
	if int(temp[index[target_point_col]]) == 0:
		continue
	chk_cnt = 0
	for neg in negates:
		if int(temp[index[neg]]) == 0:
			chk_cnt += 1
	# pass if all negates are 0(missing value)
	if chk_cnt == len(negates):
		continue

	## RULE2: All outside should be smaller (increased) or larger (decreased) that objects
	target_point_level = int(temp[index[target_point_col]])
	if tendency_direction == "increase":
		allowance_level = target_point_level - distinction_level
		chk_cnt = 0
		# pass unless all ancestors are less or eqaul to target point - distinction_level
		for anc in ancestors:
			if int(temp[index[anc]]) <= allowance_level:
				chk_cnt += 1
		if not chk_cnt == len(ancestors):
			continue

		chk_cnt = 0
		# pass unless all negates are less or equal to all ancestors
		for neg in negates:
			for anc in ancestors:
				if int(temp[index[neg]]) > int(temp[index[anc]]):
					chk_cnt += 1
		if chk_cnt > 0:
			continue
	elif tendency_direction == "decrease":
		allowance_level = target_point_level + distinction_level
		chk_cnt = 0
		# pass unless all ancestors are greater or eqaul to target point + distinction_level
		for anc in ancestors:
			if int(temp[index[anc]]) >= allowance_level:
				chk_cnt += 1
		if not chk_cnt == len(ancestors):
			continue
 
		chk_cnt = 0
		# pass unless all negates are greater or equal to all ancestors
		for neg in negates:
			for anc in ancestors:
				if int(temp[index[neg]]) < int(temp[index[anc]]):
					chk_cnt += 1
		if chk_cnt > 0:
			continue

	
	## RULE3: Sub-ancestor should be within allowed noise levels
	if not len(subancestors) == 0:
		upper = int(temp[index[target_point_col]]) + allowed_noise_leve
		lower = 0
		if (int(temp[index[target_point_col]]) - allowed_noise_level) > 0:
			lower = int(temp[index[target_point_col]]) - allowed_noise_level
		chk_cnt = 0
		for subanc in subancestors:
			if int(temp[index[subanc]]) <= upper and int(temp[index[subanc]]) >= lower:
				chk_cnt += 1
		if not chk_cnt == len(subancestors):
			continue

		""" obsolete STRICT RULE: All sub-clade has same increase/decrease tendency
		if tendency_direction == "increase":
			chk_cnt = 0
			for subanc in subancestors:
				if int(temp[index[subanc]]) < int(temp[index[target_point_col]]):
					chk_cnt += 1
			if chk_cnt > 0:
				continue
		elif tendency_direction == "decrease":
			chk_cnt = 0
			for subanc in subancestors:
				if int(temp[index[subanc]]) > int(temp[index[target_point_col]]):
					chk_cnt += 1
			if chk_cnt > 0:
				continue
		"""


	## print filtered 
	result.append(line)

f.close()


### Print out result ###
with open(outfile, 'w') as f:
	for line in result:
		f.write("%s\n" % line)
f.close()

