#! bin/bash

#Path to raw data folder
data_path="../00_raw_data/30-1072689103/00_fastq"

#Path to raw data folder containing the fastq.gz and fastq.gz.md5 files I downloaded from Genewiz.
fastq_list=$data_path/*.fastq.gz

output_path="../02_fastqc/"

if [ ! -d "$output_path" ]; then
    echo $output_path does not exist! Creating one now.
    mkdir $output_path
fi

for filename in $fastq_list; do
    fastqc $filename -o $output_path
done
