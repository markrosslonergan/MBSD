#!/bin/bash

mkdir data/angle5.0_Ratio_bound
FILENAME=data/angle5.0_Ratio_bound/bestfit.dat

echo > "$FILENAME"

#for ms in `ls ../decay/data/ | awk 'FS="_" {print $1}' | uniq | tail -n63 | sort -r -n | paste -sd ' ' -`
for ms in `seq 0.1500 -0.0100 0.0600` `seq 0.0590 -0.0020 0.0010` 
do
#if [ `echo "$ms < 0.1110" | bc -l` -eq 1 ]
#then

#for mZ in `ls ../decay/data | grep $ms"_" | awk 'FS="_" {print $2}' | sed 's/.dat//g' | uniq | sort -rn | paste -sd ' ' -`
for mZ in `seq 0.1200 0.0200 0.9800`
do
#if [ `echo "$mZ < 0.7200" | bc -l` -eq 1 ]
#then

	if [ `echo "$ms < $mZ" | bc -l` -eq 1 ];
	then 
		echo $ms $mZ
		./stats -m $ms -Z $mZ -c 5.0 -R -P >> "$FILENAME"
			
	fi
#fi
done
	echo  >> "$FILENAME"
	echo
#fi
done
