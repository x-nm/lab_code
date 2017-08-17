# Fetch genelist in source code of LifeMap webpage, data from HOMER
# Input: 
#	tissue symbol, e.g. CART-cranial-neural-crest-cells
#	url
# Output: 
#	MARKER .gmx file: MARKERS-CART-cranial-neural-crest-cells.gmx
#	EXPRESS .gmx file: EXPR-CART-cranial-neural-crest-cells.gmx
#	with url and comments in the second line, about EXPRESSION LEVEL, +/-
# 
#  4s for a request
#
# 2017.6.2 by xnm

import urllib2
import re

# INPUT
tissue = "BONE"
cell = "Membranous-Facial-Bones-Intramembranous-Preosteoblasts"
cardUrl = "in-vivo-development/bone/membranous-facial-bones/intramembranous-preosteoblasts"

symbol = tissue+"-"+cell
trunk = "https://discovery.lifemapsc.com/"
url = trunk + cardUrl

# OUTPUT
markerFile = "MARKERS-"+symbol+".gmx"
exprFile = "EXPR-"+symbol+".gmx"
outDir ="D:\\xnm\\work\\17.2.15\\KeyGene\\geneset"
file_out_marker = open(outDir+"\\"+markerFile,'w')
file_out_expr = open(outDir+"\\"+exprFile,'w')
# first line
file_out_marker.write("MARKERS-"+symbol+'\n')
file_out_expr.write("EXPR-"+symbol+"\n")

# HTTP REQUEST
response = urllib2.urlopen(url)
webpage = response.read()

# REGEX MATCH
loc4TargetText = re.compile(r'{ GeneExpression:.*}') # in one line, no need to use multi-line match,or non-greedy match...
text = loc4TargetText.findall(webpage)
text = text[0]

loc4Item = re.compile(r'{"CellID":.+?}')
#loc4Item = re.compile(r'{"TCCellID":.+?}') # FOR STEM CELL PAGES
geneSymbol = re.compile(r'(?<="GeneSymbol":")[\w,-]+')
exprLevel = re.compile(r'(?<="ExpressionLevel":)-?\d')
exprPattern = re.compile(r'(?<="ExpressionPattern":)(null|0|1|2|3)')


highSel_up = [] # highly selective genes, expression pattern = 3, level = 1
highSel_down = [] # highly selective genes, expression pattern = 3, level = -1
sel_up = [] # selective genes, expression pattern = 2, level = 1
sel_down = [] # selective genes, expression pattern = 2, level = -1
others = [] # expressed genes, not markers, expr pattern = 0/null

for item in loc4Item.findall(text):
	gene = geneSymbol.findall(item)[0].upper()
	level = exprLevel.findall(item)[0]
	pattern = exprPattern.findall(item)[0]
	if pattern == "0" or pattern == "null":
		others.append(gene)
	else:
		if pattern == "3":
			if level == "1":
				highSel_up.append(gene)
			elif level == "-1":
				highSel_down.append(gene)
		elif pattern == "2":
			if level == "1":
				sel_up.append(gene)
			else:
				sel_down.append(gene)

comment = "highDown:"+str(highSel_down)+" Down:"+str(sel_down)+"; highUp: "+str(highSel_up)
file_out_marker.write(">LifeMap Marker; "+comment+'\n')
file_out_expr.write(">LifeMap, gene expressed in the cell; "+comment+'\n')
for i in highSel_down:
	file_out_marker.write(i+"\n")
	file_out_expr.write(i+"\n")
for i in sel_down:
	file_out_marker.write(i+"\n")
	file_out_expr.write(i+"\n")
for i in highSel_up:
	file_out_marker.write(i+"\n")
	file_out_expr.write(i+"\n")
for i in sel_up:
	file_out_marker.write(i+"\n")
	file_out_expr.write(i+"\n")
for i in others:
	file_out_expr.write(i+"\n")


file_out_marker.close()
file_out_expr.close()