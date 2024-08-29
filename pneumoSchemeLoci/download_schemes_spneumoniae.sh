#!/bin/bash

wget -O loci https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/schemes/2/loci && \
cat loci | grep -o 'https[^"]*' | sed 's/$/\/alleles_fasta/' | parallel -v -k -j 5 'mkdir $(basename $(dirname {})) && wget -q -P $(basename $(dirname {})) {}'
cat SPN* >> alleles_spneumoniae.fasta && rm -rf SPN*
