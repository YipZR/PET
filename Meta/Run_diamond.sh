#/bin/bash


######################### conda source ################################

source /home/yzr/miniconda3/etc/profile.d/conda.sh

######################################################################


Path_Aspera=/home/yzr/miniconda3/etc/asperaweb_id_dsa.openssh
Path_DB=/home/yzr/Work/Meta/DB_test/ARGs_PET_MHET.dmnd
Path_SRA=/data_1/yzr_data/SRR
Path_EBI=/home/yzr/Work/Meta/EBI
Path_Kraken2=/data_1/databases/kraken2


##################

Sample='Not_ID_SRA.txt'

##################

Sample_Name=`echo $Sample | awk -F '_' '{print $1}'`
cd $Path_SRA
mkdir $Sample_Name
cd $Sample_Name
mkdir tmp
mkdir download_dir
mkdir diamond_result

###################
cd $Path_EBI
len=$( awk 'END {print NR}' $Sample )
i=2
for((;i<=$len;i++));
do
	cd $Path_EBI
        SRR=$(cat $Sample | sed -n ${i}p | awk '{print $1}')
	echo $SRR
	SRRtype=$(echo ${SRR:0:3} | tr 'A-Z' 'a-z')
	echo $SRRtype
	SRRhead=${SRR:0:6}
	SRRtail1=${SRR: -1}
	SRRtail2=${SRR: -2}
	SRRtail3=${SRR: -3}	

	url1='era-fasp@fasp.sra.ebi.ac.uk:/vol1/'$SRRtype'/'$SRRhead'/00'$SRRtail1'/'$SRR
	echo 'url1 === '$url1
	url2='era-fasp@fasp.sra.ebi.ac.uk:/vol1/'$SRRtype'/'$SRRhead'/0'$SRRtail2'/'$SRR
	echo 'url2 === '$url2
	Path_Work=$Path_SRA/$Sample_Name/
        Path_Download=$Path_Work'download_dir'
	
	cd $Path_Download
	ascp -v -QT -l 300m -P 33001 -k 1 -i $Path_Aspera $url1 ./

	NUM=$(ls ./ | wc -l)
	if ((NUM == 0)); then
		echo 'Next'        
		ascp -v -QT -l 300m -P 33001 -k 1 -i $Path_Aspera $url2 ./
	fi
        if ((NUM == 1)); then
                echo 'Download Success'
        fi
	mv ./* ..
	cd ..

	fastq-dump --split-files $SRR
	mv ./*.fastq ./tmp

	cd ./tmp
	num=$(ls ./ | grep .fastq | wc -l)
	if ((num ==1)); then
		FQ1=$(ls ./ | grep .fastq)
		FQ1Paird=(${FQ1//fastq/1Paird.fastq})
		trimmomatic SE -phred33 $FQ1 $FQ1Paird ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 -threads 100
	fi
	if ((num ==2)); then
		FQ1=$(ls ./ | grep _1.fastq)
		FQ2=$(ls ./ | grep _2.fastq)
		FQ1Paird=(${FQ1//1.fastq/1Paird.fastq})
		FQ1UnPaird=(${FQ1//1.fastq/1UnPaird.fastq})
		FQ2Paird=(${FQ1//1.fastq/2Paird.fastq})
		FQ2UnPaird=(${FQ1//1.fastq/2UnPaird.fastq})
		trimmomatic PE -phred33 $FQ1 $FQ2 $FQ1Paird $FQ1UnPaird $FQ2Paird $FQ2UnPaird ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 -threads 100
	fi
	FQclean=$(ls ./ | grep 1Paird.fastq)
#	M8=(${FQclean//1Paird.fastq/_diamond.m8})
#	M82=(${FQclean//1Paird.fastq/2_diamond.m8})
#	diamond blastx -p 144 -k 1 -d $Path_DB -q $FQclean -o $M8 
	
	conda activate kraken2
        kraken2_result=(${FQ1Paird//1Paird.fastq/_kraken2.txt})
        kraken2 --db $Path_Kraken2 $FQclean --thread 144 --use-mpa-style --report $kraken2_result

#	diamond blastx -p 100 -k 1 -d $Path_DB2 -q $FQclean -o $M82
	rm ./*.fastq
	mv ./* ..
#	conda activate mpa
#	Phlan_Out=(${FQclean//1Paird.fastq/_metaphlan.txt})
#	Bowtie_Out=(${FQclean//1Paird.fastq/.bowtie2.bz2})
	echo 'Success'
	echo $i
done
