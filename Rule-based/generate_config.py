"""
Parse newick tree and return config file.

USAGE:
python generate_config.py [newick_tree.txt] [node_name] [config_file_name]

Original: tree_parsor.py, junha.shin@wisc.edu
"""
import sys, re
from collections import defaultdict

### INPUT ###
tree_file = sys.argv[1]	# newick tree format
target = sys.argv[2]	# name of any node
outfile = sys.argv[3]	# ouput file name

f = open (tree_file, 'r')
for tree in f:
	tree = tree.rstrip()
f.close()
print (tree)
print (target)


topls = tree.split("(")
topanc = topls[0]
if not topls[0] == "Anc1":
	topls = tree.split(")")
	topanc = topls.pop()
	topls.append(topanc)
	if not topanc == "Anc1":
		sys.exit("ERROR in newick tree format: NO Anc1")
	## inverting tree
	#  step 1. split all ( , ) and nodes into array
	#  step 2. switch parenthesis direction: "( -> )", ") -> ("
	#  setp 3. invert array
	tmp1 = re.sub(r'\(','\t[\t',tree)
	tmp2 = re.sub(r',','\t,\t', tmp1)
	tmp3 = re.sub(r'\)','\t(\t',tmp2)
	tmp4 = re.sub(r'\[',')',tmp3)
	tmp5 = re.sub(r'\t+', '\t', tmp4)
	topls = tmp5.split("\t")
	topls = list(reversed(topls))
	topls.pop()
	treeinv = "".join(topls)
	#print (treeinv)
	topls = treeinv.split("(")


### Load data with depth infomation for the next step ###
#	Key idea is that since the newick tree has ordered structure, 
#	by looking at each nodes based on the order we could parse them.
alls = []		# list of all nodes
ancs = []		# list of ancestors by occurence in newick tree
anclevel = dict()	# k=ancestor_node, v=level
mat = dict()		# k=level, vARR=children ordred by occurence in newick tree
col = 0			# index for hierarchy by "(", move one step to child
for branch in topls:
	temp = branch.split(")")
	move = 0	# index for hierarchy by ")", move one step to parent
	for cell in temp:
		col = col - move
		if not cell == "":
			if mat.has_key(col):
				mat[col].append(cell)
			else:
				mat[col] = [cell]
			move = move + 1
			# no need to move when extant_end, not an ancestor-related ")"
			# therefore restore move value to before increasing
			if re.search("^,", cell):
				if not re.search("^,Anc", cell):	# means extant end
					move = move - 1
		#print (cell, col, move)
		pair = cell.split(",")
		for side in pair:
			if not side == "":
				alls.append(side)
			if re.search(r"Anc", side):
				ancs.append(side)
				anclevel[side] = col
	col = col + 1
#print(alls)
#print(ancs)
#print(anclevel)
#print(mat)


### Parse newick tree to direct children and ancestor ###
def add_parent(anc, child):
	#if direct_parent.has_key(child):
	#	direct_parent[child].append(anc)
	#else:
	#	direct_parent[child] = [anc]
	direct_parent[child] = anc

direct_children = dict()	# k=ancestor, vARR=2 children
direct_parent = dict()		# k=ancestor, v=ancestor of ancestor
for anc in ancs:
	direct_children[anc] = []
	childlevel = anclevel[anc] + 1
	childlevelarr = mat[childlevel]
	childtemp = []
	for cell in childlevelarr:
		arr = cell.split(",")
		for val in arr:
			if not val == "":
				childtemp.append(val)
	if len(childtemp) == 2:
		direct_children[anc] = childtemp
		for child in childtemp:
			add_parent(anc, child)
	else:
		# pop out first element at twice
		child = childtemp.pop(0)
		direct_children[anc].append(child)
		mat[childlevel] = childtemp
		add_parent(anc, child)
		
		child = childtemp.pop(0)
		direct_children[anc].append(child)
		mat[childlevel] = childtemp
		add_parent(anc, child)
#print (direct_children)
#print (direct_parent)


### Collect data per ancestral node ("Anc") based on parsed direct relationship ###
def collect_subanc (pair):
	for child in pair:
		if re.search(r"Anc", child):
			gathering.append(child)

def collect_ancestors (anc):
	if direct_parent.has_key(anc):
		gathering.append(direct_parent[anc])

subancs = dict()		# k=ancestor, vARR=subancestor
children = defaultdict(dict)	# k1=ancestor, k2=children, v=dummy
ancestors = defaultdict(dict)	# k1=ancestor, k2=ancestors of k1 ancestor, v=dummy
for anc in ancs:
	gathering = []
	# collect sub-ancestors
	subancs[anc] = []
	collect_subanc(direct_children[anc])
	if len(gathering) > 0:
		while len(gathering) > 0:
			subanc = gathering.pop(0)
			subancs[anc].append(subanc)
			collect_subanc(direct_children[subanc])

	# collect all children
	for child in direct_children[anc]:
		children[anc][child] = 1
	for subanc in subancs[anc]:
		children[anc][subanc] = 1
		for child in direct_children[subanc]:
			children[anc][child] = 1

	# collect all ancestors
	collect_ancestors(anc)
	if len(gathering) > 0:
		while len(gathering) > 0:
			tempanc = gathering.pop(0)
			ancestors[anc][tempanc] = 1
			collect_ancestors(tempanc)
#print (subancs)
#print (children)
#print (ancestors)


### Results ###
obj = []	# target ancestor node and its all extant incl. sub-ancestors
extant = []	# all extant (excl. ancestor nodes)
subs = []	# all sub-ancestors
ances = []	# all ancestors of target ancestor node
negate = []	# opposite of object

tgt_anc = ""    # target ancestor point
try:
	if re.search("^Anc", target):
		tgt_anc = target
	else:
		tgt_anc = direct_parent[target]
except:
	print("ERROR in target name:%s\n" % target)

## First generate sets based on Anc
obj.append(tgt_anc)
for anc, child in children.iteritems():
	if anc == tgt_anc:
		for c in child.keys():
			obj.append(c)
			if not re.search(r"Anc", c):
				extant.append(c)

for i, sc in enumerate(subancs[tgt_anc]):
	subs.append(sc)

for i, an in enumerate(ancestors[tgt_anc]):
	ances.append(an)

for node in alls:
	if not node in obj and not node in ances:
		negate.append(node)


## Deal with non-Anc target case based on firstly generated sets
if not re.search("^Anc", target):
	# update negate by adding obj except itself and direct_parent
	obj.remove(target)
	obj.remove(tgt_anc)
	for node in obj:
		negate.append(node)
	# clear and remake obj
	obj = []
	obj.append(target)
	# clear and remake extant
	extant = []
	extant.append(target)
	# just clear subs
	subs = []
	# add direct_parent into ances
	ances.append(tgt_anc)


## print out result
with open(outfile, 'w') as f:
	f.write("tree=%s\n" % tree)
	f.write("target_point=%s\n" % target)
	#objs = " ".join(obj)
	f.write("object=%s\n" % " ".join(obj))
	#exts = " ".join(extant)
	f.write("extant=%s\n" % " ".join(extant))
	if len(subs) > 0:
		f.write("subancestor=%s\n" % " ".join(subs))
	else:
		f.write("subancestor=\n")
	f.write("ancestor=%s\n" % " ".join(ances))
	f.write("negate=%s\n" % " ".join(negate))
f.close()
	
