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

    Usage: bac_ass_pipe.sh [options] FASTQ

    Options:
        -i        Input file. [required].
        -t        Threads. [4].
        -w        Working directory. Path to create the folder which will contain all mitnanex information. [./bac_ass_pipe_out].
        -r        Prefix name add to every produced file. [input file name].
        -d        Different output directory. Create a different output directory every run (it uses the date and time). [False]. 
        *         Help.
    "
    exit 1
}

while getopts 'i:t:w:r:d' opt; do
    case $opt in
        i)
        input_file=$OPTARG
        ;;
        t)
        threads=$OPTARG
        ;;
        w)
        wd=$OPTARG
        ;;
        w)
        prefix=$OPTARG
        ;;
        d)
        output_dir="mitnanex_results_$(date  "+%Y-%m-%d_%H-%M-%S")/"
        ;;
        *)
        bac_ass_pipe_help
        ;;
    esac 
done

# Check if required arguments are provided
if [ -z "$input_file" ];
then
  echo "Error: Input file is required."
  bac_ass_pipe_help
fi

## PREFIX name to use for the resulting files
if [ -z $prefix ];
then 
    prefix=$(basename $input_file)
    prefix=${prefix%%.*}
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