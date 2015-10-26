#!/bin/bash

if [ -e "/scratch/ross/dataMC12oct/anti/MC_"$1"_"$2".dat" ]; then
	echo File exists.
	./decayer $1 $2 > "data12oct/anti/"$1"_"$2".dat"
fi


