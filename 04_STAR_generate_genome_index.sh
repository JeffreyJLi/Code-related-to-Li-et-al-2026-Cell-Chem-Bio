#! bin/bash

#Path to folder containing GRCh38 human reference genome .fa and .gtf files
genome_path="../gencode_human_genome/"

#Path to GRCh38 human reference genome
genome_file_path="../gencode_human_genome/GRCh38.primary_assembly.genome.fa"

#Path to GRCh38 human reference genome .gtf annotation file
gtf_file_path="../gencode_human_genome/gencode.v46.primary_assembly.annotation.gtf"


#Path to output trimmed data.
output_path=../04_STAR_genome_indices/

#Checks if 03_trimmed output directory exists and if not, makes one.
if [ ! -d "$output_path" ]; then
    echo $output_path does not exist! Creating one now.
    mkdir $output_path
fi

# sjdbOverhang is set to the readlength - 1 (I used 2x150 bp)
STAR --runThreadN 36 --runMode genomeGenerate --genomeDir $output_path --genomeFastaFiles $genome_file_path --sjdbGTFfile $gtf_file_path --sjdbOverhang 149