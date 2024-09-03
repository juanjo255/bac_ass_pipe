#!/bin/bash
wd=$1
wget -O loci https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/schemes/2/loci && \
cat loci | grep -o 'https[^"]*' | sed 's/$/\/alleles_fasta/' | parallel -v -k -j 5 "curl --silent --create-dirs --output $wd'/scheme_alleles/'\$(basename \$(dirname {})).fa {}"
