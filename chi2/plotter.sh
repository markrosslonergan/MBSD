#!/bin/bash 

#for ms in `grep -v ^# data/bestfit_chiU.dat | awk '{print $1}' | uniq | sort -n`
#do
#for mZ in `grep ^$ms data/bestfit_chiU.dat | awk '{print $2}' | sort -rn`
#do
ms=0.0160
mZ=0.4000

	echo $ms" "$mZ
	
		CHI2_EBEST=`grep ^$ms" "$mZ" " data/bestfit_chiU.dat | awk '{print $3}'`  		
		EBEST=`grep ^$ms" "$mZ" " data/bestfit_chiU.dat | awk '{print $4}'`  		
		EBEST_CONT=`grep ^$ms" "$mZ" " data/bestfit_chiU.dat | awk '{print $5}'`  		
		CHI2_ABEST=`grep ^$ms" "$mZ" " data/bestfit_chiU.dat | awk '{print $6}'`  		
		ABEST=`grep ^$ms" "$mZ" " data/bestfit_chiU.dat | awk '{print $7}'`  		
		ABEST_CONT=`grep ^$ms" "$mZ" " data/bestfit_chiU.dat | awk '{print $8}'`  		
		CUTEFF=`grep ^$ms" "$mZ" " data/bestfit_chiU.dat | awk '{print $9}'` 
 
		MS=`printf "%.4lf" "$ms"` 
		MZ=`printf "%.4lf" "$mZ"` 
		
		./stats -m $ms -Z $mZ -cP -X $EBEST -E > "data/Espectra/Ebest_"$MS"_"$MZ".dat"
		./stats -m $ms -Z $mZ -cP -X $ABEST -A > "data/Aspectra/Abest_"$MS"_"$MZ".dat"	

		ETOT=`grep -v ^# "data/Espectra/Ebest_"$MS"_"$MZ".dat" | awk '{total = total +$2}END{print total}'`
		ATOT=`grep -v ^# "data/Aspectra/Abest_"$MS"_"$MZ".dat" | awk '{total = total +$2}END{print total}'`

		cd plots
	
		cp ballettsquartet.plt temp.plt

		FILENAME=$MS"_"$MZ

		sed -i s/FILENAME/$FILENAME/g temp.plt
		sed -i s/TEMPMS/$MS/g temp.plt
		sed -i s/TEMPMZ/$MZ/g temp.plt
		sed -i s/TEMPCUTEFF/$CUTEFF/g temp.plt
		sed -i s/TEMPECONTEFF/$EBEST_CONT/g temp.plt
		sed -i s/TEMPACONTEFF/$ABEST_CONT/g temp.plt
		sed -i s/EBESTDELTACHI/$CHI2_EBEST/g temp.plt
		sed -i s/ABESTDELTACHI/$CHI2_ABEST/g temp.plt
		sed -i s/EBESTUCHI/$EBEST/g temp.plt
		sed -i s/ABESTUCHI/$ABEST/g temp.plt
		sed -i s/TEMPETOT/$ETOT/g temp.plt
		sed -i s/TEMPATOT/$ATOT/g temp.plt



		sed -i s/EFNAME/$FILENAME/g temp.plt
		sed -i s/AFNAME/$FILENAME/g temp.plt

		gnuplot temp.plt
		epstopdf --autorotate=All $FILENAME"_quartet.eps"
		rm $FILENAME"_quartet.eps"
		rm temp.plt
		mv $FILENAME"_quartet.pdf" "pdf/"$FILENAME"_quartet.pdf"

		cd ..

#done
#	echo
#done
