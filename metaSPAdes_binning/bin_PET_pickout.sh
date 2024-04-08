#!/bin/bash

Path_out=/data_1/yzr_data/MetaSpades/OUT/bin_PET_MHET_pickout
Path_work=/data_1/yzr_data/MetaSpades/OUT


list='list.txt'

len=$( awk 'END {print NR}' $list )
echo $len
i=1
for((;i<=$len;i++));
do
        item=$(cat $list | sed -n ${i}p)
	SRR=`echo $item | awk -F '_' '{print $1}'`
	echo $SRR
	cd $SRR'_1e5'
	cd bin
	cp $item'.fa' $Path_out
	cd $Path_work
done
