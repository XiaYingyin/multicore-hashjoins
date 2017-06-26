
for ((i=1; i<=128; i=$[i*2]))
do
	echo "---------------------------------------"
	rows=$[i*600000000/100]
	percent=`printf "%0.6f\n" $(bc -q <<< scale=6\;${i}/100)`
	echo "fact rows ${rows}"
	echo "--------NPO STAR JOIN NSM--------------"
	./src/mchashjoins -a NPO_SN -n 16 -e ${percent}
	echo ""
	echo "--------NPO STAR JOIN DSM--------------"
	./src/mchashjoins -a NPO_SD -n 16 -e ${percent}
	echo ""
	echo "--------NPO STAR JOIN VECTOR--------------"
	./src/mchashjoins -a STARJOIN --starjoinflag -n 4 -e ${percent}
	echo ""
	echo "---------------------------------------"
	echo ""
	echo ""
done

