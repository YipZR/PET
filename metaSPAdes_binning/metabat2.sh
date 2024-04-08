#/bin/bash

source /home/yzr/miniconda3/etc/profile.d/conda.sh
conda activate metabat2
###################################

Path_SRRfiles=/data_1/yzr_data/SRR/landfill
Path_bin=/data_1/yzr_data/MetaSpades/OUT

####################################

Sample='MetaSpades_landfill_ID.txt'

####################################

len=$( awk 'END {print NR}' $Sample )
echo $len
i=1
for((;i<=$len;i++));
do
	SRR=$(cat $Sample | sed -n ${i}p | awk '{print $1}')
	echo $SRR
	cd $Path_SRRfiles
#	fasterq-dump -e 100 -p ./$SRR -O $Path_bin/$SRR'_1e5'
	cd $Path_bin	
	cd $SRR'_1e5'
	mkdir bin
	FQ1=$(ls ./ | grep _1.fastq)
        FQ2=$(ls ./ | grep _2.fastq)
        FQ1Paird=(${FQ1//1.fastq/1Paird.fastq})
        FQ1UnPaird=(${FQ1//1.fastq/1UnPaird.fastq})
        FQ2Paird=(${FQ1//1.fastq/2Paird.fastq})
        FQ2UnPaird=(${FQ1//1.fastq/2UnPaird.fastq})
        trimmomatic PE -phred33 $FQ1 $FQ2 $FQ1Paird $FQ1UnPaird $FQ2Paird $FQ2UnPaird ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 -threads 144

#	bowtie2-build -f $SRR'_contigs.fasta' $SRR.index --threads 100
	bowtie2 -1 *1Paird.fastq -2 *2Paird.fastq -p 144 -x $SRR.index -S $SRR.sam
	samtools view -@ 144 -b $SRR.sam -o $SRR.bam
	samtools sort -@ 144 -l 9 -O BAM $SRR.bam -o $SRR'_sorted.bam'
	jgi_summarize_bam_contig_depths --outputDepth $SRR'_depth.txt' $SRR'_sorted.bam'
	metabat2 -m 1500 -t 144 -i $SRR'_contigs.fasta' -a $SRR'_depth.txt' -o ./bin/$SRR'_bin' -v
	cd bin
	checkm lineage_wf -f checkm.txt -t 144 -x fa ./ ./checkm/ 
	grep 'bin' checkm.txt | sed 's/^  //' | awk '{print $1,$2,$13,$14}' | sed 's/\ /\t/g'| sed 's/\./\t/' | sort -n -k 2 | sed 's/\t/./' > $SRR.txt
	echo 'i ===' $i
	cd $Path_bin
done
