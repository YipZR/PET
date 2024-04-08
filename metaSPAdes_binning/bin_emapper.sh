#bin/bash

Path_blastx_out=/data_1/yzr_data/MetaSpades/OUT/blastx
Path_eggnog_out=/data_1/yzr_data/MetaSpades/OUT/eggnog
File_PETase_db=/data_1/yzr_data/MetaSpades/OUT/blastx_db/PETase.dmnd
Path_sh=/data_1/yzr_data/MetaSpades/OUT


Sample='MetaSpades_landfill_ID.txt'

len=$( awk 'END {print NR}' $Sample )
echo $len
i=14
for((;i<=$len;i++));
do
        SRR=$(cat $Sample | sed -n ${i}p | awk '{print $1}')
        echo $SRR
	cd bin_PET_pickout
        num=$(ls ./ | grep .fa | wc -l)
	j=1
	for((;j<=$num;j++));
	do
		emapper.py -i $SRR'_bin'.$j.fa --itype genome -o $Path_eggnog_out/$SRR'_bin'.$j --cpu 144 --translate -d bact
	done
	cd $Path_sh
	echo 'success'
	echo 'i === '$i
done
