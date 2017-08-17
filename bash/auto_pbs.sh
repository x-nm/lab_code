#!/bin/bash

# modify the script written in windows, 
# and auto generate a pbs.script with the simplist frame
# 2017.8.11 by xnm

USAGE="USAGE: ./auto_pbs.sh test.R"

if [ $1 == "-h" ]; then
	echo $USAGE
	exit 1
fi

# get type and filename
tp=${1##*.} # Ex. R
name=${1%.*} # Ex. test

# type 2 operation
if [ $tp == "R" ]; then
	op="Rscript"
elif [ $tp == "py" ]; then
	op="python"
else
	echo "wrong type"
	exit 1
fi

# clean the win script
cat $1 | tr -d "\r" >temp
rm $1
mv temp $1

# generate the pbs.script
cat>"$name".sh<<EOF
#!/bin/bash
#$ -cwd
$op $1
EOF


