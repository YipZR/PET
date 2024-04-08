#bin/bash
source /home/yzr/miniconda3/etc/profile.d/conda.sh

Path_sh=/data_1/yzr_data/MetaSpades/OUT
Path_bin=/data_1/yzr_data/MetaSpades/OUT/bin_PET_MHET_pickout
Path_DB=/data_1/yzr_data/MetaSpades/OUT/blastx_db


bin_list='list_PET_MHET.txt'

len=$( awk 'END {print NR}' $bin_list )
echo $len

i=1
for((;i<=$len;i++));
do
        bin=$(cat $bin_list | sed -n ${i}p | awk '{print $1}')
        echo $bin
	cd $Path_bin
	cd $bin
	seqkit grep -f $bin'_PET_ID.txt' $bin'.fa' -o $bin'_PET_contigs.fasta'
	seqkit grep -f $bin'_MHET_ID.txt' $bin'.fa' -o $bin'_MHET_contigs.fasta'
	
	conda activate prodigal
	prodigal -i $bin'_PET_contigs.fasta' -a $bin'_PET_protein.faa' -p meta
	prodigal -i $bin'_MHET_contigs.fasta' -a $bin'_MHET_protein.faa' -p meta

	diamond blastp -d $Path_DB/PETase.dmnd -q $bin'_PET_protein.faa' -f 6 -o $bin'_PET_protein.matches' --id 30 --subject-cover 80
	diamond blastp -d $Path_DB/MHETase.dmnd -q $bin'_MHET_protein.faa' -f 6 -o $bin'_MHET_protein.matches' --id 30 --subject-cover 80

	awk {'print $1'} $bin'_PET_protein.matches' | sort -u > $bin'_PET_protein_matches_ID.txt'
	awk {'print $1'} $bin'_MHET_protein.matches' | sort -u > $bin'_MHET_protein_matches_ID.txt'
	
	seqkit grep -f $bin'_PET_protein_matches_ID.txt' $bin'_PET_protein.faa' -o $bin'_final_PET.faa'
	seqkit grep -f $bin'_MHET_protein_matches_ID.txt' $bin'_MHET_protein.faa' -o $bin'_final_MHET.faa'

	cat $bin'_final_PET.faa' >> $Path_bin/binpickout_fianl_PET.faa
	cat $bin'_final_MHET.faa' >> $Path_bin/binpickout_final_MHET.faa
	
	mkdir $bin'_result'
	mv $bin'_PET_contigs.fasta' $bin'_MHET_contigs.fasta' $bin'_PET_protein.faa' $bin'_MHET_protein.faa' $bin'_PET_protein.matches' $bin'_MHET_protein.matches' $bin'_PET_protein_matches_ID.txt' $bin'_MHET_protein_matches_ID.txt' $bin'_final_PET.faa' $bin'_final_MHET.faa' $bin'_result'
	
	cd $Path_sh
	echo 'success'
	echo 'i === '$i
done
