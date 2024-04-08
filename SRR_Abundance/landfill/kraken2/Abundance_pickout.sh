

ID='landfill_Abundance_ID.txt'

Sample_Name=`echo $ID | awk -F '_' '{print $1}'`
echo $Sample_Name

len=$( awk 'END {print NR}' $ID )
echo $len

i=1
for((;i<=$len;i++));
do
        SRR=$(cat $ID | sed -n ${i}p | awk '{print $1}')
	echo $SRR

	mv $SRR'_kraken2.txt' $Sample_Name'_Abundance'
	echo $i
done
