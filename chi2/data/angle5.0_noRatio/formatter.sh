#!/bin/bash


for n in `grep -v ^# bestfit.dat | awk '{print $1}' | uniq`
do
	for m in `grep ^0.001" " bestfit.dat | awk '{print $2}'`
	do
		WORD=`grep ^$n" "$m" " bestfit.dat`

		if [ "$WORD" = "" ]
		then
			echo $n" "$m" 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
		else
			echo $WORD
		fi
	done
	echo 
done
