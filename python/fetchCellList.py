# Fetch cell list and info needed from LifeMap webpage
#
# Input: tissue name
# Output: cell info list file
# 
# 10s+ for a request
#
# 2017.06.23 by xnm

#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import urllib2

# Input
tissue = "cartilage" # bone, cartilage
trunk = "https://discovery.lifemapsc.com/in-vivo-development/"
url = trunk + tissue

# HTTP request
resonse = urllib2.urlopen(url)
webpage = resonse.read()

# REGEX match
loc4section = re.compile(r'Cells are listed according to their order of development from the first ancestors.*?],items:\[.*?\]') # because all the data are in the same line
text = loc4section.findall(webpage)
text = text[0]
# basic info
loc4Item = re.compile(r'{"ID":.+?"MatchedTCCells.*?}') # there are {} in item content...
CardUrl = re.compile(r'(?<="CardUrl":")[-,/,\w]*')
CellName = re.compile(r'(?<="Name":")[\w,\s]*')
AnatCompName = re.compile(r'(?<="AnatCompName":")[\w,\s]*')
# other info
IsPrimaryProgenitor = re.compile(r'(?<="IsPrimaryProgenitor":)\w+') # true or false
IsMainPrimaryProgenitor = re.compile(r'(?<="IsMainPrimaryProgenitor":)\w+')
TotalGeneExp = re.compile(r'(?<="TotalGeneExp":)\d*')
ValidatedGeneExp = re.compile(r'(?<="ValidatedGeneExp":)\d*')
MarkerGeneExp = re.compile(r'(?<="MarkerGeneExp":)\d*')
SignalInvolvedGeneExp = re.compile(r'(?<="SignalInvolvedGeneExp":)\d*')


# record info
outDir ="D:\\xnm\\work\\17.5.25\\LifeMap\\"
outFile = open(outDir+tissue+"_cell_info.txt",'w')
outFile.write("Name\tCardUrl\tisPrimaryProgenitor\tisMainPrimaryProgenitor\tgeneNum\tvalidatedGeneNum\tmarkerNum\tsignalInvolvedGeneNum\n")

def tf201(boolstr):
	if boolstr == "true":
		return "1"
	else:
		return "0"

for item in loc4Item.findall(text):
	card_url = CardUrl.findall(item)[0]
	Name = AnatCompName.findall(item)[0].replace(" ","-")+"-"+CellName.findall(item)[0].replace(" ","-") 
	priProgenitor = tf201(IsPrimaryProgenitor.findall(item)[0])
	mainPriProgenitor = tf201(IsMainPrimaryProgenitor.findall(item)[0])
	geneNum = TotalGeneExp.findall(item)[0]
	geneValidatedNum = ValidatedGeneExp.findall(item)[0]
	markerNum = MarkerGeneExp.findall(item)[0]
	signalGeneNum = SignalInvolvedGeneExp.findall(item)[0]
	outFile.write(Name+"\t"+card_url+"\t"+priProgenitor+"\t"+mainPriProgenitor+"\t"+geneNum+"\t"+geneValidatedNum+"\t"+markerNum+"\t"+signalGeneNum+"\n")



