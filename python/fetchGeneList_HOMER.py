# Fetch genelist in source code of LifeMap webpage, data from HOMER
# Input: tissue symbol, e.g. BONE
# Output: .gmx file: BONE.gmx
# 
# 3~7s for a request
#
# 2017.6.2 by xnm

import urllib2
import re


# INPUT
#tissue = "SKIN"
tissue = "KIDNEY"

# OUTPUT
outputFile = tissue+".gmx"
outDir ="D:\\xnm\\work\\17.2.15\\KeyGene\\geneset"
file_out = open(outDir+"\\"+outputFile,'w')
file_out.write(tissue+'\n')

# HTTP REQUEST
trunk = "https://discovery.lifemapsc.com/gene-expression-signals/high-throughput/large-scale-dataset-homer/"
#trunk = "https://discovery.lifemapsc.com/gene-expression-signals/high-throughput/rnaseq-large-scale-dataset-human-protein-atlas-rna-sequencing-data/"
url = trunk+tissue.lower()

#url = "https://discovery.lifemapsc.com/gene-expression-signals/high-throughput/large-scale-dataset-human-protein-atlas-rna-sequencing-data/kidney"

response = urllib2.urlopen(url)
webpage = response.read()

# REGEX MATCH
pattern = re.compile(r'(?<="GeneSymbol":")[\w,-]+')
match = pattern.findall(webpage)
if match:
	count = len(match)
	#file_out.write(">HOMER; "+str(count)+'\n')
	file_out.write("> "+str(count)+";"+url+'\n')
	for i in match:
		file_out.write(i+'\n')

file_out.close()