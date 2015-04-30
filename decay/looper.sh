#!/bin/bash

#echo > logfile.dat
for ms in `ls ../MC/MC_* | awk 'FS="_" {print $2}' | uniq | sort -n | paste -sd ' ' -`
do

if [ `echo "$ms-0.038 > 0" | bc -l` -eq 1 ]
then

	for mZ in `ls ../MC/MC_* | awk 'FS="_" {print $3}' | tail -n44 | sed 's/.dat//g' | sort -rn`
	do

	START=`date +%s`
	echo $ms $mZ
	./all.sh $ms $mZ
	FINISH=`date +%s`
	RUNTIME=$(( $FINISH - $START ))
	echo $ms" "$mZ" "$RUNTIME >> logfile.dat
	echo RUNTIME=$RUNTIME

	done
fi

done 

