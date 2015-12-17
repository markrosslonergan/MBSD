#!/bin/bash

if [ -e "/scratch/ross/dataMC28oct/MC_"$1"_"$2".dat" ]; then
	echo File exists.
	./decayer $1 $2 > "data28oct/"$1"_"$2".dat"
fi


