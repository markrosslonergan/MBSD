#!/bin/bash

if [ -e "/scratch/ross/latestrunSeptabove2/MC_"$1"_"$2".dat" ]; then
	echo File exists.
	./decayer $1 $2 > "data/above2/"$1"_"$2".dat"
fi


