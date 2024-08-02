#!/bin/bash

#Default values
threads=4
wd="./bac_ass_pipe_out"
memory=$(awk '/MemFree/ { printf "%.0f", $2/1024/1024 }' /proc/meminfo)

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
        -m        Memory for SKESA in GB. [$memory]
        -w        Working directory. Path to create the folder which will contain all mitnanex information. [./bac_ass_pipe_out].
        -f        FastP options. [' '].
        -n        Different output directory. Create a different output directory every run (it uses the date and time). [False].
        *         Help.
    "
    exit 1
}

while getopts '1:2:d:t:m:w:f:n' opt; do
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
        m)
        memory=$OPTARG
        ;;
        w)
        wd=$OPTARG
        echo "directory"
        ;;
        f)
        fastp_opts=$OPTARG
        ;;
        n)
        output_dir="bac_ass_pipe_results_$(date  "+%Y-%m-%d_%H-%M-%S")/"
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
    prefix1=${prefix1%%.*}
fi

if [ -z $prefix2 ];
then 
    prefix2=$(basename $R2_file)
    prefix2=${prefix2%%.*}
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
    ## Check if output directory exists
    if [ -d $1 ];
    then
        echo " "
        echo "Directory $1 exists."
        echo " "
    else 
        mkdir $1
        echo " "
        echo "Directory $1 created"
        echo " " 
    fi
}

trimming(){
    echo " " 
    echo "Step 1: Trimming using fastp"
    echo "FastP options: " $f
    echo " "
    fastp $f --thread $threads -i $R1_file -I $R2_file -o $wd$prefix1".filt.fastq.gz" -O $wd$prefix2".filt.fastq.gz" \
            -j $wd"fastp.json" -h $wd"fastp.html"
    
    # Use the filtered reads in the rest of the pipeline
        R1_file=$wd$prefix1".filt.fastq.gz"
        R2_file=$wd$prefix2".filt.fastq.gz"
}

assembly (){
    ## Unicycler
    echo "Step 2: Assemblying with Unicycler and SKESA"
    
    ## Create output folders
    create_wd $wd"skesa_asm"
    create_wd $wd"unicycler_asm"

    #unicycler -t $threads -1 $R1_file -2 $R2_file -o $wd"/unicycler_asm"
    skesa --reads $R1_file,$R2_file --cores $threads --memory $memory --contigs_out $wd"skesa_asm/assembly_skesa.fasta"
}

quality_asm (){
    echo "Step 3: Quality assessment of assembly produced using QUAST, BUSCO, Kraken2"

}

#create_wd $wd #&& trimming && assembly
#assembly

echo $wd

echo "Finished" 