#!/bin/bash

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

#Default values
threads=4
output_dir="pneumoPipe_out"
wd=$(pwd)
k=51
scaled=100
busco_dataset="lactobacillales_odb10"
path_to_scheme=$wd"/scheme_alleles_spneumoniae/"
path_to_busco_dataset=$wd
train_file_pyrodigal=$exec_path"/prodiga_training_files/prodigal_training_files/Streptococcus_pneumoniae.trn"
unicycler_asm_fasta=$unicycler_asm"/"$prefix1".fasta"

## Help message
pneumoPipe_help() {
    echo "
    pneumoPipe - Bacterial Assembly Pipeline for Illumina reads

    Author:
    Gustavo Gamez
    Juan Picon Cossio

    Version: 0.1

    Usage: 
    pneumoPipe.sh [options] -1 reads_R1.fastq -2 reads_R2.fastq
    OR
    pneumoPipe.sh [options] -3 path/to/dir/Reads

    

    Options:    
    Required:

        -1        Input R1 paired end file. [required].
        -2        Input R2 paired end file. [required].
        -3        Path to R1 and R2 files. [required].
        -d        Kraken database. if you do not have one, you need to create it first.
        -r        Reference genomes folder. This is to get the closest genome in the NCBI for a better genome quality analysis.
    
    Optional:
        -t        Threads. [4].
        -w        Working directory. Path to create the folder which will contain all pneumopipe information. [$wd].
        -o        Output directory name. If -3 option was given a new output dir name will be assigned using reads name, otherwise default. [$output_dir]
        -f        FastP options. [' '].
        -n        Different output directory. Create a different output directory every run (it uses the date and time). [False].
        -u        Update MLST and cgMLST database. Otherwise it will assume you have both databases correctly set. [False]
        -s        Path to schemes. [$path_to_scheme]
        -b        Path to busco datasets. [$path_to_busco_dataset]
        -p        Training file for pyrodigal. [$train_file_pyrodigal]. 
        *         Help.
    "
    exit 1
}

while getopts '1:2:3:d:r:t:w:o:f:nus:b:p:' opt; do
    case $opt in
        1)
        R1_file=$OPTARG
        ;;
        2)
        R2_file=$OPTARG
        ;;
        3)
        path_to_dir_paired=$OPTARG
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
        b)
        path_to_busco_dataset=$OPTARG
        ;;
         p)
        train_file_pyrodigal=$OPTARG
        ;;
        *)
        pneumoPipe_help
        ;;
    esac 
done


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
set_name_for_outfiles(){
        prefix1=$(basename $R1_file)
        prefix1=${prefix1%%.*}
        prefix2=$(basename $R2_file)
        prefix2=${prefix2%%.*}
}

if [ ${wd: -1} != / ];
then 
    wd=$wd"/"
fi

## If user decide to update database
## Skip if cgMLST db exists, otherwise download
if [ -z $update_MLST ];
then
    if [ -d $path_to_scheme ];
    then
        echo " "
        echo "Directory with cgMLST scheme exists. Skipping downloading database"
        echo " "
    else 
        mkdir "cgMLST scheme were not found at" $path_to_scheme
        echo " "
        echo "Downloading cgMLST scheme fo S. pneumoniae"
        echo " " 
        update_MLST_db
    fi
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

assembly(){
    ## Unicycler
    echo "Step 2: Assemblying with Unicycler"
    
    ## Create output folders
    unicycler_asm=$wd"unicycler_asm"
    create_wd $unicycler_asm

    unicycler --verbosity 0 -t $threads -1 $R1_file -2 $R2_file -o $unicycler_asm >> $log
    unicycler_asm_fasta=$unicycler_asm"/"$prefix1".fasta"
    mv $unicycler_asm"/assembly.fasta" $unicycler_asm_fasta

    ## Stats for report
    asm_stats=$(seqkit stats  -T --all $unicycler_asm_fasta | cut -f 4,5,6,7,8,13,17,18)

}

run_sourmash(){

    outdir_ref=$wd"/reference_genomes_signatures/"
    outdir_query=$wd"/sourmash_assess/"

    create_wd $outdir_ref
    create_wd $outdir_query

    ## Signature for query
    echo ""
    echo "Creating query signature"
    sourmash sketch dna -f -p k=$k,scaled=$scaled --outdir $outdir_query $unicycler_asm_fasta 2>> $log &&
    echo "Sourmash signature for query is at: "$outdir_query

    ## Signature for reference
    echo ""
    echo "Creating references signatures. It might take some minutes"
    sourmash sketch dna -f -p k=$k,scaled=$scaled  --outdir $outdir_ref \
        $(find $references_genomes_folder -type f -name "*.fna")  2>> $log &&
    echo "Sourmash signatures for references are at: "$outdir_ref


    ## Run search
    echo ""
    echo "Searching query signatures in reference signatures"
    sourmash search --containment -k $k $outdir_query"assembly.fasta.sig" $outdir_ref -o $outdir_query"/sourmash_out.csv" 2>> $log

}

quality_asm(){
    
    echo " "
    echo "Step 3: Quality assessment of assembly produced using QUAST, BUSCO, Kraken2 and sourmash"
    echo " "
    echo "This step will take some minutes..."      
    
    ## BUSCO UNICYCLER
    echo "Running BUSCO"
    busco -f -c $threads -m genome --download_path $path_to_busco_dataset -l $busco_dataset -i $unicycler_asm_fasta --metaeuk -o $unicycler_asm"/busco_assessment" >> $log &&

    ## SOURMASH
    echo "Running Sourmash"
    run_sourmash &&
    echo ""
    echo "Selecting the most similar reference"
    
    ## Add to report later
    reference=$(cut -d "," -f "1,3" $outdir_query"/sourmash_out.csv" | head -n 2 | grep -o "GC[^.]*") &&
    reference=$(find $references_genomes_folder -type f -name $reference"*.fna")
    reference_feature=$(grep -o ".*/" <<< $reference)"genomic.gff" 
    reference_similarity=$(cut -d "," -f "1,3" $outdir_query"/sourmash_out.csv" | head -n 2 | grep -o "^0....")

    echo "The selected reference genome is: " $reference
    echo "With a similarity of: " $reference_similarity
    
    

    ## QUAST
    echo "Running Quast"
    quast --circos --plots-format "png" -t $threads -r $reference --features "GFF:"$reference_feature -1 $R1_file -2 $R2_file \
        -o $wd"/quast_assess" $unicycler_asm_fasta 2>> $log
}

cps_serotyping(){

    echo "Step 4: Capsule serotyping using SeroCall"
    echo " "

    out_serocall=$wd"/serotype_seroCall"
    create_wd $out_serocall
    echo " "
    serocall -t $threads -o $out_serocall"/seroCall" $R1_file $R2_file &&

    ## To report
    serotype=$(cat $out_serocall"/seroCall_calls.txt" | grep -o "^[0-9].*")
}

sequence_typing (){
    
    echo " "
    echo "Step 6: Sequence typing using classic MLST and cgMLST"
    
    # create dir for MLST and cgMLST
    mlst_outdir=$wd"/mlst/"
    create_wd $mlst_outdir


    ## Classic MLST
    echo " "
    echo "Starting classic MLST"
    
    mlst --quiet --scheme "spneumoniae" --threads $threads $unicycler_asm_fasta > $mlst_outdir"MLST.tsv"
    ST=$(cat $mlst_outdir"MLST.tsv" | cut -f 3)

    ## cgMLST
    echo " "
    echo "Starting cgMLST "

    ## Prepare cgMLST scheme
    cgMLST_scheme=$wd"/cgMLST_scheme_chewBBACCA"
    chewBBACA.py PrepExternalSchema -g $path_to_scheme -o $cgMLST_scheme \
        --ptf $exec_path"/prodigal_training_files/Streptococcus_pneumoniae.trn" --cpu $threads >> $log

    ## Alellic calling
    allelic_call=$wd"/allelic_call_chewBBACCA"
    chewBBACA.py AlleleCall --cds-input -i $pyrodigal_outdir -g $cgMLST_scheme -o $allelic_call --cpu $threads --cds $pyrodigal_outdir"/genes.aa.fasta" \
        --output-novel --output-missing --no-inferred  >> $log

    ## Add to report
    ## Headers
    ### EXAC INF PLOT3 PLOT5 LOTSC NIPH NIPHEM ALM ASM PAMA LNF Invalid CDSs Classified_CDSs Total_CDSs
    cgMLST_stats=$(grep -o "\b[0-9].*" $allelic_call"/results_statistics.tsv")

}

update_MLST_db () {
    ## Update MLST and cgMLST database

    echo " "
    echo "Updating the MLST database"
    echo " "
    mlst-download_pub_mlst -j 5 -d $(echo $(grep -Po ".*(?=bin)" <<< $(which mlst))"db/pubmlst")
    mlst-make_blast_db
    
    echo "Updating the cgMLST"
    bash $exec_path"/pneumoSchemeLoci/download_schemes_spneumoniae.sh" $path_to_scheme
    echo "Files were downloaded in the working directory"
}

CDS_prediction (){
    ## CDS prediction with pyrodigal

    echo "Step 5: Gene prediction with pyrodigal"
    echo " "
    pyrodigal_outdir=$wd"/cds_pred_pyrodigal"
    create_wd $pyrodigal_outdir

    pyrodigal -j $threads -t $train_file_pyrodigal \
        -i $unicycler_asm_fasta -f "gff" -o $pyrodigal_outdir"/genes.gff" -d $pyrodigal_outdir"/genes.fasta" -p "single"
    
}

## Create report for summary of pipeline results
create_report () {

    create_wd $1 &&
    report=$wd"/report.txt"
    log=$wd"/log.txt"
    touch $report $log
    echo "-------------------------------------"
    echo "Data to proccess: " $R1_file $R2_file 

}

## Pipeline execution order
pipeline_exec(){

    #trimming && assembly && quality_asm && cps_serotyping && CDS_prediction && sequence_typing

}

export_report(){
    #ReadsID referenceID serotype_pctg_similarity ST cgMLSTStats ReferenceID referenceSimilarity assembly statistics
    echo -e "$prefix1\t$serotype\t$ST\t$cgMLST_stats\t$reference\t$reference_similarity\t$asm_stats" >> $report
    
    ## HEADERS ##
    ##headers=ReadsID referenceID serotype_pctg_similarity ST EXAC INF PLOT3 PLOT5 LOTSC \
    ## NIPH NIPHEM ALM ASM PAMA LNF Invalid CDSs Classified_CDSs Total_CDSs ReferenceID referenceSimilarity assembly statistics

}

## START PIPELINE



## Check if files or directory with files was given
## Check required files are available
if ! [ -z "$path_to_dir_paired" ];
then
    # I must attach to wd a new folder for each read
    # Here I save the common wd. I could use dirname, but this way is headache-free
    keep_wd_body=$wd
    for i in $(find $path_to_dir_paired -name "*.fastq*" | grep -o ".*R1\..*")
    do
        echo "running PneumoPipe for"  $(basename $i)
        R1_file=$i
        R2_file=$(echo "$i" | sed 's/R1/R2/')

        # Asign a name for the dir base on the reads name
        read_name=$(basename $R1_file)
        output_dir="${read_name%_*}/"

        # RUN pneumoPipe
        set_name_for_outfiles
        wd=$keep_wd_body$output_dir
        create_report $wd && 
        echo "results will be saved at" $wd
        pipeline_exec
        
    done
else
    if [ -z "$R1_file" ]; then echo "ERROR => File 1 is missing"; pneumoPipe_help; fi
    if [ -z "$R2_file" ]; then echo "ERROR => File 2 is missing"; pneumoPipe_help; fi
    
    set_name_for_outfiles
    wd=$wd$output_dir
    create_report $wd &&
    echo "results will be saved at" $wd
    pipeline_exec
fi

echo "Finished"