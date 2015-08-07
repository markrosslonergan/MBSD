#!/bin/bash

echo > logfile.dat
for ms in `ls ../MC/MC_* | awk 'FS="_" {print $2}' | uniq | sort -n | paste -sd ' ' -`
do
	for mZ in `ls ../MC/MC_"$ms"* | awk 'FS="_" {print $3}' | sed 's/.dat//g' | sort -rn`
	do

	if [ `echo "$ms-$mZ > 0" | bc -l` -eq 1 ]
	then

		START=`date +%s`
		echo $ms $mZ
		./decayer $ms $mZ > "data/2body_"$ms"_"$mZ".dat"
		FINISH=`date +%s`
		RUNTIME=$(( $FINISH - $START ))
		echo $ms" "$mZ" "$RUNTIME >> logfile.dat
		echo RUNTIME=$RUNTIME

	fi

	done
done 

