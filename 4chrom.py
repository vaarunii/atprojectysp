import sys
import gzip
import statistics
import math
import random 

# entropy calc
def entropycalc(vals):
	s = sum(vals)
	P = []
	h = 0
	for i in range(len(vals)):
		if s == 0: break
		else: p = vals[i] / s
		P.append(p)
	for p in P:
		if p == 0: continue
		h -= (p * math.log2(p))
	return h

# manhattan distance
def mandist(g1, g2):
	g1 = genex[g1]
	g2 = genex[g2]
	d = 0
	
	for tis1, tis2 in zip(g1, g2):
		dtc = abs(tis1 - tis2)
		d += dtc

	return d

# minimum difference in distance btwn two genes in cluster
def mindiff(genes):
	dtc = 2
	for i in range(len(genes)):
		for j in range(i):
			mandis = mandist(genes[j], genes[i])
			if mandis < dtc: dtc = mandis
	return dtc

genex = {}

with gzip.open(sys.argv[1], 'rt') as fp:
	for line in fp:
		gene = line.split()
		
		if gene[0] not in genex:
			if not gene[0].endswith('.1'): continue
			tissues = [0] * 11
			genex[gene[0]] = tissues

		for i in range(len(tissues)): 
			tissues[i] += int(gene[i+4])

# deleting under threshold vals
delete = []
threshold = 100
for name in genex:
	s = sum(genex[name])
	if s < threshold: delete.append(name)
for name in delete:
	del genex[name]

# converting to probabilities
for name, tissues in genex.items():
	s = sum(tissues)
	for i in range(len(tissues)):
		genex[name][i] = tissues[i] / s

# list of names
names = list(genex.keys())
names.sort()

# list of chromosome 4
chrom4 = []
for name in names:
	if name.startswith('AT4G'): chrom4.append(name)

# random pair comparison
randndn = []
for i in range(10000):
	randname1 = random.choice(chrom4)
	randname2 = random.choice(chrom4)
	if randname1 != randname2:
		expdiff = mandist(randname1, randname2)
		if expdiff not in randndn:
			randndn.append(expdiff)
	else: 
		continue
print('rand pairs:', statistics.mean(randndn), statistics.stdev(randndn))

# real next-door neighbor comparison
realndn = []
for i in range(1, len(chrom4)):
	realndn.append(mandist(chrom4[i-1], chrom4[i]))
print('real ndn:', statistics.mean(realndn), statistics.stdev(realndn))

# random neighborhoods comparison
randneighbors = []
groupsize = 20
for i in range(10000):
	randnames = []
	for j in range(groupsize):
		randname = random.choice(chrom4)
		if randname not in randnames:
			randnames.append(randname)
	randneighbors.append(mindiff(randnames))
print('rand neighbors:', statistics.mean(randneighbors), statistics.stdev(randneighbors))

# real neighborhoods comparison (slide)
realneighbors = []
groupsize = 20
for i in range(groupsize, len(chrom4)):
	realneighbors.append(mindiff(chrom4[i-groupsize:i]))
print('real neighbors (slide):', statistics.mean(realneighbors), statistics.stdev(realneighbors))

# real neighborhoods comparison (skip)
realneighbors2 = []
for i in range(groupsize, len(chrom4), groupsize):
	realneighbors2.append(mindiff(chrom4[i-groupsize:i]))
print('real neighbors (skip):', statistics.mean(realneighbors2), statistics.stdev(realneighbors2))

'''
python3 4chrom.py ~/Code/Research/DATA/at_ime_tissues.txt.gz
'''