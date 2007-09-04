#!/bin/bash

function run_slr {
	$SLRBIN -seqfile $I.paml -treefile $I.tree -outfile $I.slr_out -timemem 1 -skipsitewise 1 &> $I.slr_log
}

if [ "$1" == "" ]
then
	echo "Usage: run_failures number [binary]"
	echo "number = 0 to run all"
	exit
fi
NRUN=$1


SLRBIN=/Users/timm/src/SLR/development/slr/bin/Slr
if [ "$2" != "" ]
then
	SLRBIN=$2
fi

COUNTER=1
for I in `cat failures`
do
	echo -n "Running $I ..."
	run_slr $I
	echo " Done."
	if [ $COUNTER -eq $NRUN ]
	then
		break
	fi
	COUNTER=$[ $COUNTER + 1 ]
done
