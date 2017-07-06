
rm memory.txt
rm time.txt
echo "NPO_NSM	NPO_DSM	STAR_JOIN	PRHO_DSM	PRO_DSM" >> time.txt
echo "" >> memory.txt
for ((i=1; i<=128; i=$[i*2]))
do
	echo "---------------------------------------"
	rows=$[i*600000000/100]
	percent=`printf "%0.6f\n" $(bc -q <<< scale=6\;${i}/100)`
	
	echo "fact rows ${rows}"
	
	echo "--------NPO STAR JOIN NSM--------------"
	\time -v ./src/mchashjoins -a NPO_DSM -n 16 -e ${percent}  &> tmp.out
       	cat tmp.out
	t1=`cat tmp.out | grep "__END__" |awk 'BEGIN {m = 0} { if (m==0) m=$7; else if ($7 > m) m=$7;} END {print m}'`
	m1=`cat tmp.out | grep "Maximum resident set size (kbytes)" | awk -F ':' '{print $2}'`
	echo ""
	echo "---------------------------------------"

	echo "--------NPO STAR JOIN DSM--------------"
	\time -v ./src/mchashjoins -a NPO_DSM -n 16 -e ${percent} &> tmp.out
       	cat tmp.out
       	t2=`cat tmp.out | grep "__END__" |awk 'BEGIN {m = 0} { if (m==0) m=$7; else if ($7 > m) m=$7;} END {print m}'`
	m2=`cat tmp.out | grep "Maximum resident set size (kbytes)" | awk -F ':' '{print $2}'`
	echo ""
	echo "---------------------------------------"
	
	echo "--------NPO STAR JOIN VECTOR--------------"
	\time -v ./src/mchashjoins -a STARJOIN --starjoinflag -n 4 -e ${percent} &> tmp.out
       	cat tmp.out
       	t3=`cat tmp.out | grep "__END__" |awk 'BEGIN {m = 0} { if (m==0) m=$7; else if ($7 > m) m=$7;} END {print m}'`
	m3=`cat tmp.out | grep "Maximum resident set size (kbytes)" | awk -F ':' '{print $2}'`
	echo ""
	echo "---------------------------------------"
	
	echo "--------PRHO DSM--------------"
	\time -v ./src/mchashjoins -a PRHO_DSM -n 16 -e ${percent} &> tmp.out
       	cat tmp.out
       	t4=`cat tmp.out | grep "__END__" |awk 'BEGIN {m = 0} { if (m==0) m=$6; else if ($6 > m) m=$6;} END {print m}'`
	m4=`cat tmp.out | grep "Maximum resident set size (kbytes)" | awk -F ':' '{print $2}'`
	echo ""
	echo "---------------------------------------"
	
	echo "--------PRO DSM--------------"
	\time -v ./src/mchashjoins -a PRO_DSM -n 16 -e ${percent} &> tmp.out
       	cat tmp.out
       	t5=`cat tmp.out | grep "__END__" |awk 'BEGIN {m = 0} { if (m==0) m=$6; else if ($6 > m) m=$6;} END {print m}'`
	m5=`cat tmp.out | grep "Maximum resident set size (kbytes)" | awk -F ':' '{print $2}'`
	echo ""
	echo "---------------------------------------"
	echo "${t1}	${t2}	${t3}	${t4}	${t5}" >> time.txt
	echo "${m1}	${m2}	${m3}	${m4}	${m5}" >> memory.txt
done

