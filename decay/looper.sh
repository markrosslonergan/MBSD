#!/bin/bash

#echo > logfile.dat
for ms in `ls /scratch/ross/dataMC24sept/MC_* | awk 'FS="_" {print $2}' | uniq | sort -n | paste -sd ' ' -`
do



#if [ `echo "$ms-0.038 > 0" | bc -l` -eq 1 ]
#then
if true
then
	for mZ in `ls /scratch/ross/dataMC24sept/MC_* | uniq| awk 'FS="_" {print $3}' | sed 's/.dat//g' | uniq | sort -n| uniq `
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

