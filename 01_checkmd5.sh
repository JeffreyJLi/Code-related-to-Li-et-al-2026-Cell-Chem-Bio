#! bin/bash

###This is a basic script to cross-check the md5sum .md5 files provided by Genewiz with their associated actual .fastq.gz files.
###This is important to ensure the data was downloaded correctly.


#Path to raw data folder
data_path=/home/jli333/JL318/00_raw_data/30-1072689103/00_fastq/

#Generates list of .md5 files in raw data folder containing the fastq.gz and fastq.gz.md5 files I downloaded from Genewiz.
fastq_list=/home/jli333/JL318/00_raw_data/30-1072689103/00_fastq/*.fastq.gz.md5
#File counter to make sure all files are checked.
file_count=1

#Path to output log file that will have the results of the md5sum checks.
log_path=/home/jli333/JL318/01_md5_check/20240829_JL318_md5_check_log.txt

#Changes working directory from /scripts to /data_folder
cd $data_path

for filename in $fastq_list; do
    #Prints a terminal message to check progress and allow user to make sure all files were checked.
    echo "Checking file number $file_count"

    #Prints a header
    if (($file_count == 1)); then
        echo "Data folder checked by 01_checkmd5.sh: " $(dirname "$filename") >> $log_path
    fi

    #Runs md5sum
    md5sum -c $filename >> $log_path
    file_count=$((file_count+1))
done
