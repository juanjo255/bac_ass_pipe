#!/bin/bash

#Default values
threads=4


## Help message
bac_ass_pipe_help() {
    echo "
    bac_ass_pipe - Bacterial Assembly Pipeline for Illumina reads

    Author:
    Juan Picon Cossio

    Version: 0.1

    Usage: bac_ass_pipe.sh [options] -1 reads_R1.fastq -2 reads_R2.fastq

    Options:    
    Required:

        -1        Input R1 paired end file. [required].
        -2        Input R2 paired end file. [required].
        -d        Kraken database. if you do not have one, you need to create it first.
    
    Optional:
        -t        Threads. [4].
        -w        Working directory. Path to create the folder which will contain all mitnanex information. [./bac_ass_pipe_out].
        -f        FastP options. [" "].
        -n        Different output directory. Create a different output directory every run (it uses the date and time). [False]. 
        *         Help.
    "
    exit 1
}

while getopts '1:2:t:w:r:f:d' opt; do
    case $opt in
        1)
        R1_file=$OPTARG
        ;;
        2)
        R2_file=$OPTARG
        ;;
        d)
        kraken_db=$OPTARG
        ;;
        t)
        threads=$OPTARG
        ;;
        w)
        wd=$OPTARG
        ;;
        f)
        fastp_opts=$OPTARG
        ;;
        n)
        output_dir="mitnanex_results_$(date  "+%Y-%m-%d_%H-%M-%S")/"
        ;;
        *)
        bac_ass_pipe_help
        ;;
    esac 
done

# Check if required arguments are provided
if [ -z $R1_file ];
then
  echo "Error: Input R1 file is required."
  bac_ass_pipe_help
fi
if [ -z $R2_file ];
then
  echo "Error: Input R2 file is required."
  bac_ass_pipe_help
fi
if [ -z $kraken_db ];
then
  echo "Error: a kraken database is required."
  bac_ass_pipe_help
fi

## PREFIX name to use for the resulting files
if [ -z $prefix1 ];
then 
    prefix1=$(basename $R1_file)
    prefix1=${prefix%%.*}
fi
if [ -z $prefix2 ];
then 
    prefix2=$(basename $R2_file)
    prefix2=${prefix%%.*}
fi

if [ ${wd: -1} = / ];
then 
    wd=$wd$output_dir
else
    wd=$wd"/"$output_dir
fi

##### FUNCTIONS #####
create_wd(){
## CREATE WORKING DIRECTORY
    if [ -d $wd ]
    then
    echo $timestamp": Rewriting directory..."
    echo " "
    else 
        echo $timestamp": Creating directory..."
        echo " "
        mkdir $wd
    fi
}

trimming(){
    echo " " 
    echo "Step 1: Trimming using fastp"
    echo "FastP options: " $f
    echo " "
    fastp $f -i $R1_file -I $R2_file -o $wd$prefix1".filt.fastq.gz" -O $wd$prefix2".filt.fastq.gz"
}

trimming