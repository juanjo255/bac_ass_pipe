#!/bin/bash

#Default values
threads=4
output_dir="pneumoPipe_out"
wd=$(pwd)
k=51
scaled=100
busco_dataset="lactobacillales_odb10"
path_to_scheme=$(grep -o ".*/" <<< $wd)
path_to_busco_dataset=$(grep -o ".*/" <<< $wd)

## Help message
pneumoPipe_help() {
    echo "
    pneumoPipe - Bacterial Assembly Pipeline for Illumina reads

    Author:
    Gustavo Gamez
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
        -w        Working directory. Path to create the folder which will contain all pneumopipe information. [$wd].
        -o        Output directory name. ["pneumoPipe_out"]
        -f        FastP options. [' '].
        -n        Different output directory. Create a different output directory every run (it uses the date and time). [False].
        -u        Update MLST and cgMLST database. Otherwise it will assume you have both databases correctly set. [False]
        -s        Path to schemes. [$path_to_scheme]
        -b        Path to busco datasets. [$path_to_busco_dataset]
        *         Help.
    "
    exit 1
}

while getopts '1:2:d:r:t:w:f:nus:' opt; do
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
        w)
        wd=$OPTARG
        ;;
        f)
        fastp_opts=$OPTARG
        ;;
        n)
        output_dir="pneumoPipe_results_$(date  "+%Y-%m-%d_%H-%M-%S")/"
        ;;
        u)
        update_MLST=1
        ;;
        s)
        path_to_scheme=$OPTARG
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

    fastqc_out=$wd"QC_check"
    create_wd $fastqc_out
    fastqc --quiet --threads $threads -o $fastqc_out $R1_file $R2_file
    fastp $f --thread $threads -i $R1_file -I $R2_file -o $wd$prefix1".filt.fastq.gz" -O $wd$prefix2".filt.fastq.gz" \
            -j $wd"fastp.json" -h $wd"fastp.html"
    
    # Use the filtered reads in the rest of the pipeline
    R1_file=$wd$prefix1".filt.fastq.gz"
    R2_file=$wd$prefix2".filt.fastq.gz"

    fastqc --quiet --threads $threads -o $fastqc_out $R1_file $R2_file
    multiqc $fastqc_out -o $wd

}

assembly (){
    ## Unicycler
    echo "Step 2: Assemblying with Unicycler"
    
    ## Create output folders
    #create_wd $wd"skesa_asm"
    unicycler_asm=$wd"unicycler_asm"
    create_wd $unicycler_asm

    unicycler --verbosity 0 -t $threads -1 $R1_file -2 $R2_file -o $unicycler_asm >> $wd"/pneumoPipe.log"
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
    sourmash sketch dna -f -p k=$k,scaled=$scaled --outdir $outdir_query $unicycler_asm"/assembly.fasta" 2>> $wd"/pneumoPipe.log" &&
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
    
    echo " "
    echo "Step 3: Quality assessment of assembly produced using QUAST, BUSCO, Kraken2 and sourmash"
    echo " "
    echo "This step will take some minutes..."      
    
    ## BUSCO UNICYCLER
    echo "Running BUSCO"
    busco -f -c $threads -m genome --download_path $path_to_busco_dataset -l $busco_dataset -i $unicycler_asm"/assembly.fasta" --metaeuk -o $unicycler_asm"/busco_assessment" >> $wd"/pneumoPipe.log" &&
    ## BUSCO SKESA
    #busco -f -c $threads -m genome -l lactobacillales_odb10 -i $wd"/unicycler_asm/assembly_skesa.fasta" --metaeuk -o $wd"skesa_asm/busco_assessment"

    ## SOURMASH
    echo "Running Sourmash"
    run_sourmash &&
    echo ""
    echo "Selecting the most similar reference"
    reference=$(cut -d "," -f "1,3" $outdir_query"/sourmash_out.csv" | head -n 2 | grep -o "GC[^.]*") &&
    reference=$(find $references_genomes_folder -type f -name $reference"*.fna")
    reference_feature=$(grep -o ".*/" <<< $reference)"genomic.gff" 
    reference_similarity=$(cut -d "," -f "1,3" $outdir_query"/sourmash_out.csv" | head -n 2 | grep -o "^0....")
    echo "The selected reference genome is: " $reference
    echo "With a similarity of: " $reference_similarity
    
    ## Add to report
    echo "The selected reference genome is: " $reference >> $report
    echo "With a similarity of: " $reference_similarity >> $report

    ## QUAST
    echo "Running Quast"
    quast --circos --plots-format "png" -t $threads -r $reference --features "GFF:"$reference_feature -1 $R1_file -2 $R2_file \
        -o $wd"/quast_assess" $unicycler_asm"/assembly.fasta" 2>> $wd"/pneumoPipe.log"
}

cps_serotyping (){

    echo "Step 4: Capsule serotyping using SeroCall"
    echo " "

    out_serocall=$wd"/serotype_seroCall"
    create_wd $out_serocall
    echo " "
    serocall -t $threads -o $out_serocall"/seroCall" $R1_file $R2_file &&

    ## Add to report
    echo "cps serotyping" >> $report
    cat $out_serocall"/seroCall_calls.txt" >> $report
}

sequence_typing (){
    
    echo " "
    echo "Sequence typing using classic MLST and cgMLST"

    ## If user decide to update database
    if [ -z $update_MLST ];
    then
    update_MLST_db
    fi

    ## Classic MLST
    echo " "
    echo "Starting classic MLST"
    mlst --quiet --scheme "spneumoniae" --threads $threads $unicycler_asm"/assembly.fasta" > $wd"MLST.tsv"
    cat $wd"MLST.tsv" >> $report

    ## cgMLST
    echo " "
    echo "Starting cgMLST "
    cgMLST_scheme=$wd"/cgMLST_scheme_chewBBACCA"
    create_wd $cgMLST_scheme
    
    ## Prepare cgMLST scheme
    chewBBACA.py PrepExternalSchema -g $path_to_scheme -o $cgMLST_scheme --ptf $exec_path"/prodigal_training_files/Streptococcus_pneumoniae.trn" --cpu $threads >> $report

    allelic_call=$wd"/allelic_call"
    ## Alellic calling
    chewBBACA.py AlleleCall -i $unicycler_asm -g $cgMLST_scheme -o $allelic_call --cpu $threads --output-novel --output-missing --no-inferred >> $report

}

update_MLST_db () {

    echo " "
    echo "Updating the MLST database"
    echo " "
    mlst-download_pub_mlst -j 5 -d $(echo $(grep -o ".*(?=bin)" <<< $(which mlst))"db/pubmlst")
    mlst-make_blast_db
    
    echo "Updating the cgMLST"
    bash $exec_path"/pneumoSchemeLoci/download_schemes_spneumoniae.sh" $path_to_scheme
    echo "Files were downloaded one directory before the working directory"
}

## Create report for summary of pipeline results
create_wd $wd &&
report=$wd"report.txt"
touch $report
echo "-------------------------" >> $report
echo "Data to proccess: " $R1_file $R2_file  >> $report


## START PIPELINE


## Executable path
# FIXME 
# This part could generate problems, I am trying just to get the executable path to locate other folders
if which pneumoPipe.sh > /dev/null 2>&1; then
    ## This wont work for a while I guess. This is thought for a wet dream of using in anaconda or nextflow
    exec_path=$(grep -o ".*/" $(which pneumoPipe.sh))
else
    if [ "$(grep -o "/" <<< $0 | wc -l)" -gt 0 ]; then
        exec_path=$(grep -o ".*/" <<< $0)
    else
        exec_path="./"
    fi
fi

#trimming && assembly && quality_asm && cps_serotyping && sequence_typing

update_MLST_db

echo "Finished"