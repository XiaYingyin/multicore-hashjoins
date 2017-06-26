echo "------  6000000 rows"
./src/mchashjoins -a NPO_SN -n 8 -e 0.01
./src/mchashjoins -a NPO_SD -n 8 -e 0.01
./src/mchashjoins -a STARJOIN --starjoinflag -n 2 -e 0.01
echo "--------------------"

echo "------  60000000 rows"
./src/mchashjoins -a NPO_SN -n 8 -e 0.10
./src/mchashjoins -a NPO_SD -n 8 -e 0.10
./src/mchashjoins -a STARJOIN --starjoinflag -n 2 -e 0.05
echo "--------------------"

echo "------  600000000 rows"
./src/mchashjoins -a NPO_SN -n 8 -e 1.0
./src/mchashjoins -a NPO_SD -n 8 -e 1.0
./src/mchashjoins -a STARJOIN --starjoinflag -n 2 -e 1.0
echo "--------------------"



