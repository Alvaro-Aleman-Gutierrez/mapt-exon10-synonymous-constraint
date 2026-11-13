#!/usr/bin/env bash
# Paso 10a: descarga genomas de Hyracoidea y crea bases BLAST de nucleótidos.
set -euo pipefail
cd "$(dirname "$0")/../data/genomes"
for acc in GCA_004026925.3 GCA_028571685.1; do
  [ -s "$acc.fna" ] || {
    curl -sSL -o "$acc.zip" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA,PROT_FASTA,GENOME_GFF"
    unzip -o -q "$acc.zip" -d "$acc"
    cat "$acc"/ncbi_dataset/data/"$acc"/*.fna > "$acc.fna"
  }
  [ -s "$acc.nsq" ] || [ -s "$acc.00.nsq" ] || [ -s "$acc.nal" ] || makeblastdb -in "$acc.fna" -dbtype nucl -parse_seqids -out "$acc" > "$acc.makeblastdb.log"
done
ls -la
