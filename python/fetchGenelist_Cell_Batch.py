# Fetch genelist in source code of LifeMap webpage: 
# 	cell name and url are input from collected info, a file
#
# Input:
#	a file contains: 
#		tissue symbol and cell name, e.g. CART, cranial-neural-crest-cells
#		url
#		Name	CardUrl	isPrimaryProgenitor	isMainprimaryProgenitor	geneNum	validatedGeneNum	markerNum	signalInvolvedGeneNum	Tissue
#		1,2;3,4;9 (-1)
# Output: 
#	MARKER .gmx file: MARKERS-CART-cranial-neural-crest-cells.gmx
#	EXPRESS .gmx file: EXPR-CART-cranial-neural-crest-cells.gmx
#	with url and comments in the second line, about EXPRESSION LEVEL, +/-
# 	A LIST OF CELL-NAMES processed
#
#  ABOUT 30s for a request
#
# 2017.6.2 by xnm
# 2017.7.4 modified for batch

import urllib2
import re

trunk = "https://discovery.lifemapsc.com/"


def fetchGenelist(tissue, cell, cardUrl):
	symbol = tissue+"-"+cell
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
	file_out_expr.write(">"+url+";LifeMap, gene expressed in the cell; "+comment+'\n')
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

# INPUT
filedir = "D:\\xnm\\work\\17.5.25\\LifeMap\\"
cellFile = filedir+"Bone_Cart_CellInfo.txt"
cellList = open(cellFile,'r')
cellList.readline() #skip the first line

# OUTPUT
cellPrcs = open(filedir+"cell_processed.txt",'w')
for line in cellList:
	line = line.strip()
	data = line.split("	")
	tis = data[8].upper() # tissue
	name = data[0] # cell
	curl = data[1] # cardUrl
	isPP = data[2]
	isMPP = data[3]
	if isPP == '1' and isMPP == '1':
		title = tis+"-"+name
		cellPrcs.write("EXPR-"+title+"	"+trunk+curl+'\n')
		cellPrcs.write("MARKERS-"+title+"	"+trunk+curl+'\n')
		fetchGenelist(tis, name, curl)
		
cellList.close()
cellPrcs.close()

