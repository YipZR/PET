#bin/bash

Path_blastx_out=/data_1/yzr_data/MetaSpades/OUT/blastx_MHETase
Path_eggnog_out=/data_1/yzr_data/MetaSpades/OUT/eggnog
File_MHETase_db=/data_1/yzr_data/MetaSpades/OUT/blastx_db/DB_Orthosearch/DB_MHETase.dmnd
Path_sh=/data_1/yzr_data/MetaSpades/OUT


Sample='MetaSpades_landfill_ID.txt'

len=$( awk 'END {print NR}' $Sample )
echo $len
i=1
for((;i<=$len;i++));
do
        SRR=$(cat $Sample | sed -n ${i}p | awk '{print $1}')
        echo $SRR
	cd $SRR'_1e5'
	cd bin
        num=$(ls ./ | grep .fa | wc -l)
	j=1
	for((;j<=$num;j++));
	do
		diamond blastx -d $File_MHETase_db -q $SRR'_bin.'$j.fa -o $Path_blastx_out/$SRR'_bin.'$j.matches -f 6 -p 50 --subject-cover 80
	done
	cd $Path_sh
	echo 'success'
	echo 'i === '$i
done
