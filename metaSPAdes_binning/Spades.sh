#!/bin/bash
################## conda source ####################

source /home/yzr/miniconda3/etc/profile.d/conda.sh

####################################################

################### Path Setting #####################

Path_SRA=/data_1/yzr_data/MetaSpades/landfill
Path_ID=/home/yzr/Work/MetaSpades
Path_HMM=/home/yzr/Work/Ortho/result/result/
Path_out=/data_1/yzr_data/MetaSpades/OUT/
################## Input Your ID file #####################

Sample='MetaSpades_landfill_ID.txt'
HMM1='IsPETase.hmm'
HMM2='IsMHETase.hmm'

###########################################################


#Sample_Name=`echo $Sample | awk -F '_' '{print $1}'`
HMM_Name1=`echo $HMM1 | awk -F '.' '{print $1}'`
HMM_Name2=`echo $HMM2 | awk -F '.' '{print $1}'`
echo $HMM_Name1
echo $HMM_Name2

len=$( awk 'END {print NR}' $Sample )


i=1
for((;i<=$len;i++));
do
	cd $Path_ID
        SRR=$(cat $Sample | sed -n ${i}p | awk '{print $1}')
	echo $SRR

	cd $Path_SRA
        mkdir $SRR'_out'
	fasterq-dump -e 100 -p ./$SRR -O ./

	conda activate spades
	metaspades.py -t 120 -k 21,33,55,77,99,127 --pe1-1 *_1.fastq --pe1-2 *_2.fastq -o $SRR'_out'
	rm ./*fastq

	cd ./$SRR'_out'
	mkdir $SRR'_hmm'
	mkdir $SRR'_prodigal'

	conda activate prodigal
	prodigal -i contigs.fasta -a protein_seq.fasta -d nucleotide_seq.fasta -o genes.gff -s poteintial.stat
	cp protein_seq.fasta ./$SRR'_hmm'
	mv protein_seq.fasta nucleotide_seq.fasta genes.gff poteintial.stat ./$SRR'_prodigal'

	cd ./$SRR'_hmm'
	cp $Path_HMM/$HMM1 ./
	cp $Path_HMM/$HMM2 ./	
	hmmsearch $HMM1 protein_seq.fasta > $SRR'_'$HMM_Name1'hmm.out'
	hmmsearch $HMM2 protein_seq.fasta > $SRR'_'$HMM_Name2'hmm.out'
	
	cd ..
	mv contigs.fasta $SRR'_contigs.fasta'
	mv $SRR'_contigs.fasta' $SRR'_hmm' $SRR'_prodigal' $Path_out

	cd $Path_SRA
	rm -rf $SRR'_out'

	echo 'Success'
	echo 'i ==== '$i
done	
