#!/bin/bash

fasta_file=$1
phmmer_database_file=$2
hmmsearch_database_file=$3

phmmer -o /dev/null -A phmmer_results.sto $fasta_file $phmmer_database_file
hmmbuild -o /dev/null --amino phmmer_results.hmm phmmer_results.sto
hmmsearch -o /dev/null -A hmmsearch_results.sto --domtblout hmmsearch_results.out phmmer_results.hmm $hmmsearch_database_file

python run-fasta.py phmmer_results.hmm
