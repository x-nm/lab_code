##!/bin/bash
# cannot use in 202.120.45.100
############################################################
# Rabbit project: from raw counts to DEG
#
# Input: 
#	a file listing labels of groups for comparation
# Output: 
#	for each comparation group: 
#		DEG_ALL_INFO: up.csv, down.csv
#		DEG_LIST: up.txt, down.txt; <gene symbol, gene name>
#		MA-plot;
#	DEG_count.txt
# Usage: RabbitDEG.sh op_file |tee report.txt
#
# Scripts used:
#	RabbitDEG.sh
#	RabbitDeseq2.R
#	Ann.py
#
# Input folder:
#	raw data file (4 groups); (defined in the R script)
#	scripts used; (defined in this script)
#	gene symbol annotation file; (defined in this script)
#	gene name annotation file; (defined in this script)
#	README.txt: for naming conventions
#
# Writer: xnm
# Created at 2016.09.09
############################################################

USAGE="USAGE: ./$0 <file for operation instructions> |tee report.txt"

if [ $# -ne 1 ]; then
        echo $USAGE
        exit 1
fi

# Data&Scripts preparation
op=$1 # the 1st arg is the operation file
DEG_R=RabbitDeseq2.R
x_ann=Ann.py
ann_symbol=gene_id2symbol.txt
ann_name=gene_symbol2name.txt

if [ ! -d output ]; then
	mkdir output
fi

# main process
IFS=" "
while read G1 G2; do
	# 0. statement
	label=DEG_$G1$G2
	echo "################################"
	echo \>""$label"":
	echo "################################"
	# 1. deseq2
	echo "STEP1: DESeq2"
	Rscript $DEG_R $G1 $G2
	# generates two files: DEG_$G1$G2_up/down.txt;
	# 2. annotation: gene symbol
	echo "STEP2: GENE SYMBOL ANNOTATION"
	python $x_ann "$label"_up.txt $ann_symbol 1 sy
	python $x_ann "$label"_down.txt $ann_symbol 1 sy
	# generates two files: DEG_$G1$G2_up/down_sy.txt;
	# 3. annotation: gene name
	echo "STEP3: GENE NAME ANNOTATION"
	python $x_ann "$label"_up_sy.txt $ann_name 8 nm
	python $x_ann "$label"_down_sy.txt $ann_name 8 nm
	# generates two files: DEG_$G1$G2_up/down_sy_nm.txt;
	# 4. cut -> list
	echo "STEP4: CUT TO LIST"
	cut -f 8,9 "$label"_up_sy_nm.txt >DL_"$G1$G2"_up.txt
	cut -f 8,9 "$label"_down_sy_nm.txt >DL_"$G1$G2"_down.txt
	# generates two files: DL_$G1$G2_up/down.txt, 
	# with format "gene symbol TAB gene name"
	# 5. DEG counting
	echo "STEP5: DEG counting"
	# usage of function: count up
	DEG_count(){
		ct=`wc -l "$label"_$1.txt`
		ct=${ct% *} #forward cut-out, as in some wc there will be filename after count
		echo "$label"_$1: $((ct-1))
	}
	DEG_count up >> DEG_count.txt
	DEG_count down >> DEG_count.txt
	# 6. redundancy removing & move into output folder
	echo "STEP6: ENDING"
	rm "$label"_up.txt
	rm "$label"_down.txt
	rm "$label"_up_sy.txt
	rm "$label"_down_sy.txt
	mv "$label"_up_sy_nm.txt "$label"_up.txt
	mv "$label"_down_sy_nm.txt "$label"_down.txt
	# left 4 txt files: 2 DEG with ann; 2 DL(DEG list)
	mv *"$G1$G2"* output 
done < $op

