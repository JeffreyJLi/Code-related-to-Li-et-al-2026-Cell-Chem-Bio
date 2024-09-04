#! bin/bash

#Path to raw data folder
data_path=/home/jli333/JL318/00_raw_data/30-1072689103/00_fastq/

#Path to raw data folder containing the fastq.gz and fastq.gz.md5 files I downloaded from Genewiz.
fastq_list=/home/jli333/JL318/00_raw_data/30-1072689103/00_fastq/*.fastq.gz
#File counter to make sure all files are checked.

output_path=/home/jli333/JL318/02_fastqc/

if [ ! -d "$output_path" ]; then
    echo $output_path does not exist! Creating one now.
    mkdir $output_path
fi

for filename in $fastq_list; do
    fastqc $filename -o $output_path
done
