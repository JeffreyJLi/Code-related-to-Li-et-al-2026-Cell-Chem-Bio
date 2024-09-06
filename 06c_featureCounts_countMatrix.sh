#! bin/bash

# This script uses featureCounts to count the number of reads per gene. However, this script passes ALL the .bam files into featureCounts at the same time, generating a count matrix that is more easily parsed by DESeq2.

#Path to raw data folder
star_path=../05_STAR_alignments

#Path to Gencode .gtf file
gtf_path="../gencode_human_genome/gencode.v46.primary_assembly.annotation.gtf"

#Declare array variable containing all 24 library labels.
declare -a condition_list=( "JL318-Rep1-A" "JL318-Rep1-B" "JL318-Rep1-C" "JL318-Rep1-D" "JL318-Rep1-E" "JL318-Rep1-F" "JL318-Rep1-G" "JL318-Rep1-H" "JL318-Rep2-A" "JL318-Rep2-B" "JL318-Rep2-C" "JL318-Rep2-D" "JL318-Rep2-E" "JL318-Rep2-F" "JL318-Rep2-G" "JL318-Rep2-H" "JL318-Rep3-A" "JL318-Rep3-B" "JL318-Rep3-C" "JL318-Rep3-D" "JL318-Rep3-E" "JL318-Rep3-F" "JL318-Rep3-G" "JL318-Rep3-H" )

#Path to output data.
output_path=../06c_featureCounts_countMatrix

#Counter to keep track of progress.
counter=1

#Checks if output directory exists and if not, makes one.
if [ ! -d "$output_path" ]; then
    echo $output_path does not exist! Creating one now.
    mkdir $output_path
fi

bam_file_list=$star_path/*.sortedByCoord.out.bam

output_filename=$output_path/JL318_featureCounts_countMatrix.txt

#Iterates through sample condition list as declared in the condition_list variable
featureCounts -T 48 -p --countReadPairs -s 1 -Q 10 -a $gtf_path -o $output_filename $bam_file_list
