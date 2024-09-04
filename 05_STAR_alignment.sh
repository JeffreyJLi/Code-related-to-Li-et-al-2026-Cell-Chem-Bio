#! bin/bash

#Path to raw data folder
trimmed_path=../03_trimmed/

#Path to raw data folder containing the fastq.gz and fastq.gz.md5 files I downloaded from Genewiz.
fastq_list=../03_trimmed/*.trimmed.fastq.gz

#Path to STAR genome indices directory
genome_index_path="../04_STAR_genome_indices/"

#Declare array variable containing all 24 library labels.
declare -a condition_list=( "JL318-Rep1-A" "JL318-Rep1-B" "JL318-Rep1-C" "JL318-Rep1-D" "JL318-Rep1-E" "JL318-Rep1-F" "JL318-Rep1-G" "JL318-Rep1-H" "JL318-Rep2-A" "JL318-Rep2-B" "JL318-Rep2-C" "JL318-Rep2-D" "JL318-Rep2-E" "JL318-Rep2-F" "JL318-Rep2-G" "JL318-Rep2-H" "JL318-Rep3-A" "JL318-Rep3-B" "JL318-Rep3-C" "JL318-Rep3-D" "JL318-Rep3-E" "JL318-Rep3-F" "JL318-Rep3-G" "JL318-Rep3-H" )

#Path to output trimmed data.
output_path=../05_STAR_alignments/

#Path to log directory
log_path=$output_path"log/"

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
    echo "Running STAR on paired end files."
    #Generates full file names of Read1 and Read2 fastq.gz files
    read1_file=$trimmed_path$condition"_R1_001.trimmed.fastq.gz"
    read2_file=$trimmed_path$condition"_R2_001.trimmed.fastq.gz"

    #Generates the full input argument
    star_input_arg=$read1_file" "$read2_file

    #Generates output file filename
    output_filename=$output_path$condition"_STAR"

    #Generates log filename
    log_filename=$log_path$condition"_log.txt"

    if [ -f $read1_file ] && [ -f $read2_file ]; then
        #Argument explanations:
        #runThreadN = number of cores to be used
        #runMode = what mode should STAR use (I want to align reads)
        #genomeDir = directory where the STAR index (which I built in 04_STAR_generate_genome_index.sh) is
        #readFilesIn = where are the .fastq.gz read1 and read2 files to be aligned?
        #readFilesCommand = required since the input files are compressed as .gz files
        #outFilterMultimapNmax = how many times can a read map to be counted as a valid alignment (lower is more stringent)
        #outFilterType = only keep the reads that contain junctions that passed filtering to SJ.out.tab
        #outSAMattributes = list which SAM attributes (we want all)
        #outSAMtype = what output format do you want? We want BAM sorted by coordinate
        #outFileNamePrefix = output files prefix
        #quantMode = Count reads per gene while mapping

        STAR --runThreadN 36 --runMode alignReads --genomeDir $genome_index_path --readFilesIn $star_input_arg --readFilesCommand zcat --outFilterMultimapNmax 1 --outFilterType BySJout --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $output_filename --quantMode GeneCounts > $log_filename
    else
        echo "ERROR: "$condition".fastq.gz paired-end files not found!"
    fi
    counter=$((counter+1))
    echo "-------------------"
done