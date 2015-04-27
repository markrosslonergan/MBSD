#!/bin/bash

echo > logfile.dat
for ms in `ls ../MC/MC_* | awk 'FS="_" {print $2}' | uniq | sort -n | paste -sd ' ' -`
do

for mZ in `ls ../MC | grep "MC_"$ms | awk 'FS="_" {print $3}' | sed 's/.dat//g' | uniq | sort -rn | paste -sd ' ' -`
do

	START=`date +%s`
	echo $ms $mZ
	./all.sh $ms $mZ
	FINISH=`date +%s`
	RUNTIME=$(( $FINISH - $START ))
	echo $ms" "$mZ" "$RUNTIME >> logfile.dat
	echo RUNTIME=$RUNTIME

done

done 

