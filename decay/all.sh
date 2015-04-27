#!/bin/bash

if [ -e "../MC/MC_"$1"_"$2".dat" ]; then
	echo File exists.
	./decayer $1 $2 > "data/"$1"_"$2".dat"
fi


