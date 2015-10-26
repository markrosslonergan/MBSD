#!/bin/bash

mkdir data/12oct
#FILENAME=data/afterdelete/parallelE1.dat

#for ms in `ls ../decay/data/ | awk -F'_' ' {print $1}' | uniq | sort -n `
#do
#parallel -j4 ./stats -m $ms -Z {1} -c 20.0 -R -P -S 0.03 -IE >> "$FILENAME" ::: `ls ../decay/data | grep $ms* | while read line; do awk -F "_" '{print $2}' <<< "$line"; done | sed 's/.dat//g' | sort -n | uniq`

#echo  >> "$FILENAME"
#done

#FILENAME=data/24sept/pizero_parallelA.dat

#for ms in `ls ../decay/data/ | awk -F'_' ' {print $1}' | uniq | sort -n `
#do
#parallel -j4  ./stats -m $ms -Z {1} -c 20.0 -R -P -S 0.03 -IA >> "$FILENAME" ::: `ls ../decay/data | grep $ms* | while read line; do awk -F "_" '{print $2}' <<< "$line"; done | sed 's/.dat//g' | sort -n | uniq`

#echo  >> "$FILENAME"
#done


FILENAME=data/12oct/gofE_dual_a.dat

for ms in `ls ../decay/data12oct/anti/ | awk -F'_' ' {print $1}' | uniq | sort -n `
do
parallel -j4  ./stats -m $ms -Z {1} -c 20.0 -R -P -S 0.01 -IE >> "$FILENAME" ::: `ls ../decay/data12oct/anti/ | grep $ms* | while read line; do awk -F "_" '{print $2}' <<< "$line"; done | sed 's/.dat//g' | sort -n | uniq`

echo  >> "$FILENAME"
done


