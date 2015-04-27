#!/bin/bash

FILENAME=data/angle5.0_noRatio/bestfit.dat

echo > "$FILENAME"
echo > logfile.dat
for ms in `ls ../decay/data/ | awk 'FS="_" {print $1}' | uniq | tail -n63 | sort -r -n | paste -sd ' ' -`
do
#ms=0.0210

#if [ `echo "$ms < 0.1110" | bc -l` -eq 1 ]
#then

for mZ in `ls ../decay/data | grep $ms"_" | awk 'FS="_" {print $2}' | sed 's/.dat//g' | uniq | sort -rn | paste -sd ' ' -`
do
	if [ `echo "$ms < $mZ" | bc -l` -eq 1 ];
	then 
		echo $ms $mZ >> logfile.dat
		echo $ms $mZ
		./stats -m $ms -Z $mZ -c 5.0 -P >> "$FILENAME"
#		./stats $ms $mZ 0.0 180.0
			
	fi
done
	echo  >> "$FILENAME"
	echo
#fi
done
