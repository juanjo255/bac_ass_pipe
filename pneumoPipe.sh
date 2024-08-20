#!/bin/bash

#Default values
threads=4
wd="./pneumoPipe_out"
memory=$(awk '/MemFree/ { printf "%.0f", $2/1024/1024 }' /proc/meminfo)
k=51
scaled=100

## Help message
pneumoPipe_help() {
    echo "
    pneumoPipe - Bacterial Assembly Pipeline for Illumina reads

    Author:
    Juan Picon Cossio

    Version: 0.1

    Usage: pneumoPipe.sh [options] -1 reads_R1.fastq -2 reads_R2.fastq

    Options:    
    Required:

        -1        Input R1 paired end file. [required].
        -2        Input R2 paired end file. [required].
        -d        Kraken database. if you do not have one, you need to create it first.
        -r        Reference genomes folder. This is to get the closest genome in the NCBI for a better genome quality analysis.
    
    Optional:
        -t        Threads. [4].
        -m        Memory for SKESA in GB. [$memory]
        -w        Working directory. Path to create the folder which will contain all mitnanex information. [./pneumoPipe_out].
        -f        FastP options. [' '].
        -n        Different output directory. Create a different output directory every run (it uses the date and time). [False].
        *         Help.
    "
    exit 1
}

while getopts '1:2:d:r:t:m:w:f:n' opt; do
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
        r)
        references_genomes_folder=$OPTARG
        ;;
        t)
        threads=$OPTARG
        ;;
        m)
        memory=$OPTARG
        ;;
        w)
        wd=$OPTARG
        ;;
        f)
        fastp_opts=$OPTARG
        ;;
        n)
        output_dir="pneumoPipe_results_$(date  "+%Y-%m-%d_%H-%M-%S")/"
        ;;
        *)
        pneumoPipe_help
        ;;
    esac 
done

# Check if required arguments are provided
if [ -z $R1_file ];
then
  echo "Error: Input R1 file is required."
  pneumoPipe_help
fi

if [ -z $R2_file ];
then
  echo "Error: Input R2 file is required."
  pneumoPipe_help
fi

# if [ -z $kraken_db ];
# then
#   echo "Error: a kraken database is required."
#   pneumoPipe_help
# fi

if [ -z $references_genomes_folder ];
then
  echo "Error: a folder with reference genomes is required."
  pneumoPipe_help
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
    echo "Step 2: Assemblying with Unicycler"
    
    ## Create output folders
    #create_wd $wd"skesa_asm"
    create_wd $wd"unicycler_asm"

    unicycler -t $threads -1 $R1_file -2 $R2_file -o $wd"/unicycler_asm" > $wd"/pneumoPipe.log"
    #skesa --reads $R1_file,$R2_file --cores $threads --memory $memory --contigs_out $wd"skesa_asm/assembly_skesa.fasta"
}

run_sourmash(){

    outdir_ref=$wd"/reference_genomes_signatures/"
    outdir_query=$wd"/sourmash_assess/"

    create_wd $outdir_ref
    create_wd $outdir_query

    ## Signature for query
    echo ""
    echo "Creating query signature"
    sourmash sketch dna -f -p k=$k,scaled=$scaled --outdir $outdir_query $wd"/unicycler_asm/assembly.fasta" 2>> $wd"/pneumoPipe.log" &&
    echo "Sourmash signature for query is at: "$outdir_query

    ## Signature for reference
    echo ""
    echo "Creating references signatures. It might take some minutes"
    sourmash sketch dna -f -p k=$k,scaled=$scaled  --outdir $outdir_ref \
        $(find $references_genomes_folder -type f -name "*.fna")  2>> $wd"/pneumoPipe.log" && 
    echo "Sourmash signatures for references are at: "$outdir_ref


    ## Run search
    echo ""
    echo "Searching query signatures in reference signatures"
    sourmash search --containment -k $k $outdir_query"assembly.fasta.sig" $outdir_ref -o $outdir_query"/sourmash_out.csv" 2>> $wd"/pneumoPipe.log"

}

quality_asm (){
    
    echo ""
    echo "Step 3: Quality assessment of assembly produced using QUAST, BUSCO, Kraken2 and sourmash"
    
    ## BUSCO UNICYCLER
    busco -f -c $threads -m genome -l lactobacillales_odb10 -i $wd"/unicycler_asm/assembly.fasta" --metaeuk -o $wd"/unicycler_asm/busco_assessment" 2>> $wd"/pneumoPipe.log" &&
    ## BUSCO SKESA
    #busco -f -c $threads -m genome -l lactobacillales_odb10 -i $wd"/unicycler_asm/assembly_skesa.fasta" --metaeuk -o $wd"skesa_asm/busco_assessment"

    ## SOURMASH
    run_sourmash &&
    echo ""
    echo "Selecting the most similar reference"
    reference=$(cut -d , -f 1,3 $outdir_query"/sourmash_out.csv" | head -n 2 | grep -o "GC[^.]*") &&
    reference=$(find $references_genomes_folder -type f -name $reference"*.fna")
    echo "The selected reference genome is: " $reference

    ## QUAST
    quast -t $threads -r $reference -1 $R1_file -2 $R2_file -o $wd"/quast_assess" $wd"/unicycler_asm/assembly.fasta" 2>> $wd"/pneumoPipe.log"
}

cps_serotyping (){

    echo "Step 4: Capsule serotyping using SeroCall"
    echo ""
    out_serocall=$wd"/serotype_seroCall"
    create_wd $out_serocall

    echo ""
    serocall -t $threads -o $out_serocall"/seroCall" $R1_file $R2_file
}

create_wd $wd && trimming && assembly && quality_asm && cps_serotyping

echo "Finished"