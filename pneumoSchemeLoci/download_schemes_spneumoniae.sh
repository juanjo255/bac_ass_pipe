#!/bin/bash
wd=$1
wget -O $wd"/loci" https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/schemes/2/loci && \
cat $wd"/loci" | grep -o 'https[^"]*' | sed 's/$/\/alleles_fasta/' | parallel -v -k -j 5 "curl --silent --create-dirs --output $wd'/scheme_alleles_spneumoniae/'\$(basename \$(dirname {})).fa {}"
