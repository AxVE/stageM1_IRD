#!/bin/bash

# Author: Tom Harrops
# Edit	: Axel Verdier

############      SGE CONFIGURATION      ###################
# Join error
#$ -j y

# Which shell
#$ -S /bin/bash 

# Export environment variables
#$ -V

# Email pour suivre l'execution 
#$ -M axel.verdier@ird.fr

# Type de massage que l'on reçoit par mail
#    -  (b) un message au demarrage
#    -  (e) a la fin
#    -  (a)  en cas d'abandon
#$ -m bea 

# Queue que l'on veut utiliser
#$ -q bioinfo.q

# Nom du job
#$ -N subsamp

############################################################


path_to_tmp="/scratch/averdier-subsample-"$JOB_ID""

file1="nas2:/data/projects/evoreprice/Data/OJ/Trimmed/OJ-cutadapt_cluster-20150522-9firstBasesTrimmed/J1_AGTCAA_L002_R1.fastq.fastq.gz"
file2="nas2:/data/projects/evoreprice/Data/OJ/Trimmed/OJ-cutadapt_cluster-20150522-9firstBasesTrimmed/J1_AGTCAA_L002_R2.fastq.fastq.gz"

outputDir="nas2:/data/projects/evoreprice/Data_test/subsample/J1"

nbReads_in_file="null" #0 : automatique
nbReads_out=1000000

# Parameters parsing
while [ $# -ge 1 ]
do
	case $1 in
		-1|--file1)
			file1=$2
			shift
		;;
		-2|--file2)
			file2=$2
			shift
		;;
		-n|--nbReadsInFile)
			nbReads_in_file=$2
			shift
		;;
		-o|--outDir)
			outputDir=$2
			shift
		;;
	esac
	shift
done





# Tmp dir creation and transfert
echo -e "\n[ $(date) : files download ]"
mkdir -p "$path_to_tmp"
cd $path_to_tmp
scp $file1 $path_to_tmp/R1.fastq.gz
scp $file2 $path_to_tmp/R2.fastq.gz
echo -e "[ $(date): done ]"

# Count nb reads if user didn't give it
if [ "$nbReads_in_file" == "null" ]; then
	echo -e "\n[$(date): count reads]"
	nbReads_in_file=$(fastq-stats R1.fastq.gz | awk '/reads\t[0-9]+/ {print $2}') #Get the number of reads by parsing the file with fastq-stats
	echo -e "[ $(date): done ]"
fi

# Display some informations
echo -e "\n==========\n[ Informations ]" 
echo -e "\tLaunch :\t$(date)"
echo -e "\tName : \t${JOB_NAME}"
echo -e "\tNode :\t${HOSTNAME}"
echo -e "\tID :\t${JOB_ID}"
echo -e "\tUser:\t${LOGNAME}"

frac=$(echo "$nbReads_out/$nbReads_in_file" | bc -l)

echo -e "\n\tnbReadsIn :\t${nbReads_in_file}"
echo -e "\tnbReadsOut :\t${nbReads_out}"
echo -e "\tfrac : \t${frac}"
echo -e "==========\n"

# Unzip files
echo -e "\n[ $(date): gunzip fastq files ]"

gunzip -c R1.fastq.gz > R1.fastq
gunzip -c R2.fastq.gz > R2.fastq

echo -e "[ $(date): done ]"
ls -lh "$path_to_tmp"

echo -e "\n[ $(date): subsampling ]"
/usr/local/bin/python /data/projects/evoreprice/Scripts/subsampleCreation/subsampler_PE_anders.py \
	$frac "$path_to_tmp"/R1.fastq "$path_to_tmp"/R2.fastq\
	$path_to_tmp/R1_test.fastq\
	$path_to_tmp/R2_test.fastq
echo -e "[ $(date): done ]"

echo -e "\n[ $(date): upload files ]"
	mkdir -p $outputDir
	scp -rp $path_to_tmp/R1_test.fastq $path_to_tmp/R2_test.fastq $outputDir

echo -e "\n[ $(date): done, cleaning up ]"
cd ../
rm -r "$path_to_tmp"
exit 0

