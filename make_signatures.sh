#!/bin/bash

k="51"
scaled="100"
common_path="/home/labcompjavier/juanPicon/neumococos_gustavo/references_genomes/asm_comparison_sourmash/"

## Signature for query
sourmash sketch dna -f -p k=$k,scaled=$scaled  --outdir $common_path "/home/labcompjavier/juanPicon/neumococos_gustavo/bac_ass_pipe_out/skesa_asm/assembly_skesa.fasta"

sourmash sketch dna -f -p k=$k,scaled=$scaled  --outdir $common_path"/references_signatures/" \
	$(find "/home/labcompjavier/juanPicon/neumococos_gustavo/references_genomes/genomes/" -type f -name "*.fna")
