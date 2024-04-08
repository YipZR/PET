#bin/bash

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
	cd *hmm
	awk '{print $9}' $SRR'_IsMHETasehmm_1e2.out' > $SRR'_IsMHETasehmm_1e2_ID.txt'
	seqkit grep -f $SRR'_IsMHETasehmm_1e2_ID.txt' protein_seq.fasta -o $SRR'_IsMHETasehmm_1e2_seq.fasta'
	seqkit seq $SRR'_IsMHETasehmm_1e2_seq.fasta' -n -o $SRR'_IsMHETasehmm_1e2_seqID.txt'
	grep 'partial=00' $SRR'_IsMHETasehmm_1e2_seqID.txt' > $SRR'_IsMHETasehmm_1e2_00ID.txt'
	seqkit grep -f $SRR'_IsMHETasehmm_1e2_00ID.txt' -n protein_seq.fasta -o $SRR'_IsMHETasehmm_1e2_00seq.fasta'
	cat $SRR'_IsMHETasehmm_1e2_00seq.fasta'>>$Path_sh/MHET_merged_00seq.fasta
	cd $Path_sh
	echo 'success'
	echo 'i === '$i
done
