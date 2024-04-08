#!/bin/bash
PATH_Query=/home/yzr/Work/Ortho/query
PATH_Genome=/home/data/yzr_data/Genomes
PATH_OrthoResult=/home/yzr/Work/Ortho/OrthoResult
len=$( awk 'END {print NR}' faa_files_bacteria.txt)

for ((i=1;i<=$len;i++));
do
	URL=$(sed -n ${i}p faa_files_bacteria.txt)
	wget -c $URL -P ./query
	cd ./query/
	if [ -e *.gz ]
	then
		gunzip *.gz
	        mv ./*.faa $i.faa	
		orthofinder -f ./
		cd ./*/*/Single_Copy_Orthologue_Sequences
		for file in ./*
		do
			if [ "${file##*.}"x = "fa"x ]
			then
				filename=`basename $file`
				temp_filename=`basename $file  .fa`
				suf=Summarize_$i
				new_filename=${suf}
				echo $filename >> ${new_filename}
				cmd="awk '{if(/>/ || /M/)print}' ${filename} >> ${new_filename}"
				eval $cmd
				echo '====== Success ======='
			else
				echo '====== 0 Result ======'
				rm -rf $PATH_Query/OrthoFinder/
			fi
		done
	else
		echo "==================================================="
		echo "=================No Protein.faa===================="
		echo "==================================================="
	fi
	mv Summarize_$i /home/yzr/Work/Ortho/Summary/$i
        mv $PATH_Query/OrthoFinder/ $PATH_OrthoResult/$i
	cd $PATH_Query
	mv ./*.faa $PATH_Genome
	cd ..
	echo 'i ==== '$i
done
