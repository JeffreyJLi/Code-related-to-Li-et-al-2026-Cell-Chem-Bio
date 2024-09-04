#! bin/bash

# Path to STAR alignments folder
star_path=../05_STAR_alignments/

#Declare array variable containing all 24 library labels.
declare -a condition_list=( "JL318-Rep1-A" "JL318-Rep1-B" "JL318-Rep1-C" "JL318-Rep1-D" "JL318-Rep1-E" "JL318-Rep1-F" "JL318-Rep1-G" "JL318-Rep1-H" "JL318-Rep2-A" "JL318-Rep2-B" "JL318-Rep2-C" "JL318-Rep2-D" "JL318-Rep2-E" "JL318-Rep2-F" "JL318-Rep2-G" "JL318-Rep2-H" "JL318-Rep3-A" "JL318-Rep3-B" "JL318-Rep3-C" "JL318-Rep3-D" "JL318-Rep3-E" "JL318-Rep3-F" "JL318-Rep3-G" "JL318-Rep3-H" )

#Counter to keep track of progress.
counter=1

#Iterates through sample condition list as declared in the condition_list variable
for condition in "${condition_list[@]}"; do
    echo "Condition "$counter": "$condition
    echo "Tabulating uniquely mapped reads for $condition"
    #Generates full file names of STAR log files
    log_file=$star_path$condition"_STARLog.final.out"

    if [ -f $log_file ]; then
        python tabulate_uniquely_mapped_reads_05b.py $log_file
    else
        echo "ERROR: "$condition" STARlog.final.out file not found!"
    fi
    counter=$((counter+1))
    echo "-------------------"
done