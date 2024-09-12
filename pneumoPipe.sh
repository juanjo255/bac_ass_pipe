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
path_to_GPSC_db=$wd"/GPSC_db/"
link_to_GPSC_db="https://gps-project.cog.sanger.ac.uk/GPS_v9.tar.gz"
link_to_GPSC_meta="https://gps-project.cog.sanger.ac.uk/GPS_v9_external_clusters.csv"
train_file_pyrodigal=$exec_path"/prodiga_training_files/prodigal_training_files/Streptococcus_pneumoniae.trn"
unicycler_outDir_fasta=$unicycler_outDir"/"$prefix1".fasta"
ariba_db="card"
ariba_db_out="card.db"
qfile=$wd"/qfile.txt"

## Help message
pneumoPipe_help() {
    echo "
    pneumoPipe - Bacterial Assembly Pipeline for Illumina reads

    Authors:
    Juan Picón Cossio
    Gustavo Gámez

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

## Add last slash
if [ ${wd: -1} != / ];
then 
    wd=$wd"/"
fi

## outDIrs
outdirs(){
    fastqc_outDir=$wd"QC_check"
    unicycler_outDir=$wd"unicycler_outDir"
    sourmash_outDir_ref=$wd"/reference_genomes_signatures/"
    sourmash_outDir_query=$wd"/sourmash_assess/"
    serocall_outDir=$wd"/serotype_seroCall"
    mlst_outDir=$wd"/mlst/"
    allelic_call_outDir=$wd"/allelic_call_outDir_chewBBACCA"
    pyrodigal_outDir=$wd"/cds_pred_pyrodigal"
    ariba_outDir=$wd"/ariba_out/"
    poppunk_outDir=$wd"/poppunk_out/"
}

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


##### FUNCTIONS #####
aesthetics(){
    echo ""
    printf '#%.0s' {1..100}
    echo ""
}

create_wd(){
    ## CREATE WORKING DIRECTORY
    ## Check if output directory exists
    if [ -d $1 ];
    then
        echo " "
        echo "Directory $1 exists."
        echo " "
    else 
        mkdir -p $1
        echo " "
        echo "Directory $1 created"
        echo " " 
    fi
}

trimming(){
    aesthetics
    echo "Step 1: Trimming using fastp"
    echo "FastP options: " $f
    aesthetics

    create_wd $fastqc_outDir
    fastqc --quiet --threads $threads -o $fastqc_outDir $R1_file $R2_file
    fastp $f --thread $threads -i $R1_file -I $R2_file -o $wd$prefix1".filt.fastq.gz" -O $wd$prefix2".filt.fastq.gz" \
            -j $wd"fastp.json" -h $wd"fastp.html"
    
    # Use the filtered reads in the rest of the pipeline
    R1_file=$wd$prefix1".filt.fastq.gz"
    R2_file=$wd$prefix2".filt.fastq.gz"

    fastqc --quiet --threads $threads -o $fastqc_outDir $R1_file $R2_file
    multiqc $fastqc_outDir -o $wd

}

assembly(){
    ## Unicycler
    aesthetics
    echo "Step 2: Assemblying with Unicycler"
    aesthetics

    ## Create output folders
    create_wd $unicycler_outDir

    unicycler --verbosity 0 -t $threads -1 $R1_file -2 $R2_file -o $unicycler_outDir >> $log
    unicycler_outDir_fasta=$unicycler_outDir"/"$prefix1".fasta"
    mv $unicycler_outDir"/assembly.fasta" $unicycler_outDir_fasta

    ## Stats for report
    asm_stats=$(seqkit stats  -T --all $unicycler_outDir_fasta | cut -f 4,5,6,7,8,13,17,18)

}

run_sourmash(){

    create_wd $sourmash_outDir_ref
    create_wd $sourmash_outDir_query

    ## Signature for query
    echo ""
    echo "Creating query signature"
    sourmash sketch dna -f -p k=$k,scaled=$scaled --outdir $sourmash_outDir_query $unicycler_outDir_fasta 2>> $log &&
    echo "Sourmash signature for query is at: "$sourmash_outDir_query

    ## Signature for reference
    echo ""
    echo "Creating references signatures. It might take some minutes"
    sourmash sketch dna -f -p k=$k,scaled=$scaled  --outdir $sourmash_outDir_ref \
        $(find $references_genomes_folder -type f -name "*.fna")  2>> $log &&
    echo "Sourmash signatures for references are at: "$sourmash_outDir_ref


    ## Run search
    echo ""
    echo "Searching query signatures in reference signatures"
    sourmash search --containment -k $k $sourmash_outDir_query"assembly.fasta.sig" $sourmash_outDir_ref -o $sourmash_outDir_query"/sourmash_out.csv" 2>> $log

}

quality_asm(){
    
    aesthetics
    echo "Step 3: Quality assessment of assembly produced using QUAST, BUSCO, Kraken2 and sourmash"
    aesthetics
    echo "This step will take some minutes..."      
    
    ## BUSCO UNICYCLER
    echo "Running BUSCO"
    busco -f -c $threads -m genome --download_path $path_to_busco_dataset -l $busco_dataset -i $unicycler_outDir_fasta -o $unicycler_outDir"/busco_assessment" >> $log &&

    ## SOURMASH
    echo "Running Sourmash"
    run_sourmash &&
    echo ""
    echo "Selecting the most similar reference"
    
    ## Add to report later
    reference=$(cut -d "," -f "1,3" $sourmash_outDir_query"/sourmash_out.csv" | head -n 2 | grep -o "GC[^.]*") &&
    reference=$(find $references_genomes_folder -type f -name $reference"*.fna")
    reference_feature=$(grep -o ".*/" <<< $reference)"genomic.gff" 
    reference_similarity=$(cut -d "," -f "1,3" $sourmash_outDir_query"/sourmash_out.csv" | head -n 2 | grep -o "^0....")

    echo "The selected reference genome is: " $reference
    echo "With a similarity of: " $reference_similarity
    
    

    ## QUAST
    echo "Running Quast"
    quast --circos --plots-format "png" -t $threads -r $reference --features "GFF:"$reference_feature -1 $R1_file -2 $R2_file \
        -o $wd"/quast_assess" $unicycler_outDir_fasta 2>> $log
}

cps_serotyping(){

    aesthetics
    echo "Step 4: Capsule serotyping using SeroCall"
    aesthetics

    create_wd $serocall_outDir
    echo " "
    serocall -t $threads -o $serocall_outDir"/seroCall" $R1_file $R2_file &&

    ## To report
    serotype=$(cat $serocall_outDir"/seroCall_calls.txt" | grep -o "^[0-9].*")
}

sequence_typing (){
    
    aesthetics
    echo "Step 6: Sequence typing using classic MLST and cgMLST"
    aesthetics

    # create dir for MLST and cgMLST
    mlst_outDir=$wd"/mlst/"
    create_wd $mlst_outDir


    ## Classic MLST
    echo " "
    echo "Starting classic MLST"
    
    mlst --quiet --scheme "spneumoniae" --threads $threads $unicycler_outDir_fasta > $mlst_outDir"MLST.tsv"
    ST=$(cat $mlst_outDir"MLST.tsv" | cut -f 3)

    ## cgMLST
    echo " "
    echo "Starting cgMLST "

    ## Prepare cgMLST scheme
    cgMLST_scheme=$wd"/cgMLST_scheme_chewBBACCA"
    chewBBACA.py PrepExternalSchema -g $path_to_scheme -o $cgMLST_scheme \
        --ptf $exec_path"/prodigal_training_files/Streptococcus_pneumoniae.trn" --cpu $threads >> $log

    ## Alellic calling
    chewBBACA.py AlleleCall --cds-input -i $pyrodigal_outDir -g $cgMLST_scheme -o $allelic_call_outDir --cpu $threads --cds $pyrodigal_outDir"/genes.aa.fasta" \
        --output-novel --output-missing --no-inferred  >> $log

    ## Add to report
    ## Headers
    ### EXAC INF PLOT3 PLOT5 LOTSC NIPH NIPHEM ALM ASM PAMA LNF Invalid CDSs Classified_CDSs Total_CDSs
    cgMLST_stats=$(grep -o "\b[0-9].*" $allelic_call_outDir"/results_statistics.tsv")

}

CDS_prediction (){
    ## CDS prediction with pyrodigal
    
    aesthetics
    echo "Step 5: Gene prediction with pyrodigal"
    aesthetics
    
    create_wd $pyrodigal_outDir

    pyrodigal -j $threads -t $train_file_pyrodigal \
        -i $unicycler_outDir_fasta -f "gff" -o $pyrodigal_outDir"/genes.gff" -a $pyrodigal_outDir"/genes.aa.fasta" -d $pyrodigal_outDir"/genes.fasta" -p "single"
    
    ## Headers are wrong in the protein file of prodigal
    echo "Fixing protein headers"
    sed 's/#.*//g; s/>\([^ ]*\)/& ID=\1;/' $pyrodigal_outDir"/genes.aa.fasta" > $pyrodigal_outDir"/genes_fixed.aa.fasta"

}

AMR_and_virulence(){

    aesthetics
    echo "Step 7: AMR and virulence factor annotation using ARIBA and AMRFinderPlus"
    aesthetics

    create_wd $ariba_outDir
    if [ -f $ariba_outDir$ariba_db_out".fa" && -f $ariba_outDir$ariba_db_out".tsv" ];
    then
        echo ""
        echo "Database for ARIBA exists, skipping download"
        echo ""
    else
        echo ""
        echo "Downloading " $ariba_db "database"
        echo ""

        ariba getref $ariba_outDir$ariba_db $ariba_outDir$ariba_db_out
        ariba prepareref --threads $threads -f $ariba_outDir$ariba_db_out".fa" -m $ariba_outDir$ariba_db_out".tsv" $ariba_outDir"ariba.$ariba_db.prepareref"
    fi

    echo ""
    echo "Running ARIBA pipeline"
    echo ""
    ariba run $ariba_outDir"ariba.$ariba_db.prepareref" $R1_file $R2_file $ariba_outDir"ariba.run"

}

GPSC_assign(){
    create_wd $poppunk_outDir
    poppunk_assign --threads $threads --db $path_to_GPSC_db"GPS_v9" --external-clustering $path_to_GPSC_db"GPS_v9_external_clusters.csv" \
        --query $qfile --output $poppunk_outDir --update-db
    #poppunk_visualise --grapetree --ref-db $path_to_GPSC_db"GPS_v9" --threads $threads --output $poppunk_outDir"/viz_all_db"
    #poppunk_visualise --grapetree --ref-db $path_to_GPSC_db"GPS_v9" --query-db --threads $threads --output $poppunk_outDir"/viz_all_db"

}

## Create report for summary of pipeline results
create_report () {

    report=$wd"/report.txt"
    log=$wd"/log.txt"
    touch $report $log

}

## Pipeline execution order for each sample
pipeline_exec_per_strain(){

    trimming && assembly && quality_asm && cps_serotyping && CDS_prediction && sequence_typing

}

## Pipeline execution for all samples previuosly proccessed

pipeline_exec_all_strain(){
    GPSC_assign

}
export_to_report(){
    #SampleID referenceID serotype_pctg_similarity ST cgMLSTStats ReferenceID referenceSimilarity assembly statistics
    echo -e "$sampleID\t$serotype\t$ST\t$cgMLST_stats\t$reference\t$reference_similarity\t$asm_stats" >> $report
    
    ## HEADERS ##
    ##headers=sampleID referenceID serotype_pctg_similarity ST EXAC INF PLOT3 PLOT5 LOTSC \
    ## NIPH NIPHEM ALM ASM PAMA LNF Invalid CDSs Classified_CDSs Total_CDSs ReferenceID referenceSimilarity assembly statistics

    ## This part creates the qfile.txt that popPUNK for GPSC assignemnt needs
    echo -e "$sampleID\t$unicycler_outDir_fasta" >> $qfile
}

## START PIPELINE

## Check if files or directory with files was given
## Check required files are available
if ! [ -z "$path_to_dir_paired" ];
then
    # I must attach to wd a new folder for each read
    # Here I save the common wd. I could use dirname, but this way is headache-free
    keep_wd_body=$wd

    ## Create one report for all samples
    create_report
    for i in $(find $path_to_dir_paired -name "*.fastq*" | grep -o ".*R1\..*")
    do
        R1_file=$i
        R2_file=$(echo "$i" | sed 's/R1/R2/')

        sampleID=$(echo "$i" | sed 's/R1.*//')

        aesthetics
        echo "Running pneumoPipe for: " $sampleID
        aesthetics
        # Asign a name for the dir base on the reads name
        read_name=$(basename $R1_file)
        output_dir="${read_name%_*}/"

        # RUN pneumoPipe
        set_name_for_outfiles
        wd=$keep_wd_body$output_dir
        create_wd $wd
        ## Create outdirs with new wd
        outdirs
        echo "results will be saved at" $wd
        pipeline_exec_per_strain
        export_to_report
    done
else
    if [ -z "$R1_file" ]; then echo "ERROR => File 1 is missing"; pneumoPipe_help; fi
    if [ -z "$R2_file" ]; then echo "ERROR => File 2 is missing"; pneumoPipe_help; fi
    
    aesthetics
    echo "Running pneumoPipe for: " $R1_file $R2_file 
    aesthetics
    set_name_for_outfiles
    wd=$wd$output_dir
    create_report
    create_wd $wd
    ## Create outdirs with new wd
    outdirs
    echo "results will be saved at" $wd
    pipeline_exec_per_strain
    export_to_report
fi
echo "Finished"