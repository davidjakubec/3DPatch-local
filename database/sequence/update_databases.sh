#!/bin/bash

for file in uniprot_sprot.fasta pdb_seqres_prot.txt; do
    if [ -e $file ]; then
        rm $file
    fi
done

wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/pdb/derived_data/pdb_seqres.txt.gz
gunzip pdb_seqres.txt.gz

grep "mol:protein" -A 1 pdb_seqres.txt > pdb_seqres_prot.txt

rm pdb_seqres.txt
