#!/bin/bash

### initial download ###

rsync -rlptvz --delete rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/mmCIF .

for mmCIF in $(find ./mmCIF/* | grep ".cif.gz"); do
    python3 generate_structure_residue_scheme.py $mmCIF
done

### update ###
