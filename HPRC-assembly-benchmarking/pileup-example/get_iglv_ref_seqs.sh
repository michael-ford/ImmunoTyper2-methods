#! /usr/bin/env bash

DB="../../ImmunoTyper2/immunotyper/data/allele_databases/IGLV/IGLV-IMGT-allele-db.fa"

for allele in 01 03 04 05; do
  grep -A1 -F "IGLV2-14*${allele}" "$DB" \
    | awk '/^>/{ print; next }{ print toupper($0) }' \
    > "iglv2-14_${allele}.fasta"
done
