#!/bin/bash

## Update MLST, cgMLST and GPSC db
wd=$1
path_to_scheme=$wd"/scheme_alleles_spneumoniae/"
path_to_GPSC_db=$wd"/GPSC_db/"
link_to_GPSC_db="https://gps-project.cog.sanger.ac.uk/GPS_v9.tar.gz"
link_to_GPSC_meta="https://gps-project.cog.sanger.ac.uk/GPS_v9_external_clusters.csv"
ariba_db="card"
ariba_db_out="card.db"
ariba_outDir=$wd"/ariba_out/"
 

echo "Updating the MLST database"
mlst-download_pub_mlst -j 5 -d $(echo $(grep -Po ".*(?=bin)" <<< $(which mlst))"db/pubmlst")
mlst-make_blast_db

echo "Downloading the cgMLST"
bash $exec_path"/pneumoSchemeLoci/download_schemes_spneumoniae.sh" $path_to_scheme

echo "Downloading ARIBA db"
mkdir -p $ariba_outDir
ariba getref $ariba_outDir$ariba_db $ariba_outDir$ariba_db_out
ariba prepareref -f $ariba_outDir$ariba_db_out".fa" -m $ariba_outDir$ariba_db_out".tsv" $ariba_outDir"ariba.$ariba_db.prepareref"

echo "Downloading GPSC db and metadata (10GB)"
mkdir -p $path_to_GPSC_db
wget -O $path_to_GPSC_db"$(basename $link_to_GPSC_db)" $link_to_GPSC_db
wget -O $path_to_GPSC_db"$(basename $link_to_GPSC_meta)" $link_to_GPSC_meta



echo "Files were downloaded in $wd"