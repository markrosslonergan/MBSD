#!/bin/bash

if [ -e "/scratch/ross/dataMC24sept/MC_"$1"_"$2".dat" ]; then
	echo File exists.
	./decayer $1 $2 > "data/"$1"_"$2".dat"
fi


