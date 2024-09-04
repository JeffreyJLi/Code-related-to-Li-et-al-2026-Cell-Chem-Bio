#! bin/bash

#Path to raw data folder
data_path=../00_raw_data/30-1072689103/00_fastq

#Path to raw data folder containing the fastq.gz and fastq.gz.md5 files I downloaded from Genewiz.
fastq_list=$data_path/*.fastq.gz

#Declare array variable containing all 24 library labels.
declare -a condition_list=( "JL318-Rep1-A" "JL318-Rep1-B" "JL318-Rep1-C" "JL318-Rep1-D" "JL318-Rep1-E" "JL318-Rep1-F" "JL318-Rep1-G" "JL318-Rep1-H" "JL318-Rep2-A" "JL318-Rep2-B" "JL318-Rep2-C" "JL318-Rep2-D" "JL318-Rep2-E" "JL318-Rep2-F" "JL318-Rep2-G" "JL318-Rep2-H" "JL318-Rep3-A" "JL318-Rep3-B" "JL318-Rep3-C" "JL318-Rep3-D" "JL318-Rep3-E" "JL318-Rep3-F" "JL318-Rep3-G" "JL318-Rep3-H" )

#Path to output trimmed data.
output_path=../03_trimmed/

#Path to output cutadapt log directory
log_path=../03_trimmed/logs/

#Counter to keep track of progress.
counter=1

#Checks if 03_trimmed output directory exists and if not, makes one.
if [ ! -d "$output_path" ]; then
    echo $output_path does not exist! Creating one now.
    mkdir $output_path
fi

#Checks if 03_trimmed/logs/ output directory exists and if not, makes one.
if [ ! -d "$log_path" ]; then
    echo $log_path does not exist! Creating one now.
    mkdir $log_path
fi

#Iterates through sample condition list as declared in the condition_list variable
for condition in "${condition_list[@]}"; do
    echo "Condition "$counter": "$condition
    echo "Running Cutadapt on paired end files."
    #Generates full file names of Read1 and Read2 fastq.gz files
    read1_file=$data_path$condition"_R1_001.fastq.gz"
    read2_file=$data_path$condition"_R2_001.fastq.gz"

    #Generates full output names of trimmed Read1 and Read 2 fastq.gz files
    trimmed_read1_file=$output_path$condition"_R1_001.trimmed.fastq.gz"
    trimmed_read2_file=$output_path$condition"_R2_001.trimmed.fastq.gz"

    #Generates summary file name
    summary_file=$log_path$condition"_cutadapt_summary.txt"

    if [ -f $read1_file ] && [ -f $read2_file ]; then
        #NOTE: Had to add -u 3 (--cut) argument to remove the first three nucleotides of the first sequencing read (Read 1), which is derived from the template-switching oligo as described in the TakaraBio library prep manual.
        #NOTE: I also had to add --minimum-length=20 since I was initially running into issues with STAR not allowing reads of zero length. 
        cutadapt --cores=48 -u 3 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length=20 -o $trimmed_read1_file -p $trimmed_read2_file $read1_file $read2_file > $summary_file
    else
        echo "ERROR: "$condition".fastq.gz paired-end files not found!"
    fi
    counter=$((counter+1))
    echo "-------------------"
done