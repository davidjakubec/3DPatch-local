# 3DPatch-tools

## About

TBD

## Requirements

    - Python >=3.5 (for `subprocess.run()`)
    - [Requests](http://docs.python-requests.org/en/master/#) Python 3 library (for making [Skylign](http://skylign.org/) requests)
    - [HMMER v3.1b2](http://hmmer.org/) program binaries in `$PATH`
    - `phmmer`/`hmmsearch` target sequence database files; examples can be downloaded using `database/sequence/update_databases.sh`
    - Up-to-date PDB structure residue scheme files; can be generated using `database/PDB_mmCIFs/update_mmCIFs.sh`
    - `run_fasta_local.py` requires `run_hmm_local.py` be present in the same directory
    - Write permission in the working directory

## Using 3DPatch-tools

### Profile HMM input

Use the file `code/run_hmm_local.py` as

    python3 run_hmm_local.py PROFILE_HMM_FILE HMMSEARCH_TARGET_DATABASE_FILE STRUCTURE_RESIDUE_SCHEMES_DIRECTORY

The `HMMSEARCH_TARGET_DATABASE_FILE` should be a file containing sequences of PDB chains. Assuming this file and the PDB structure residue schemes have been generated as described in the [Requirements](#requirements), the following example should work

    mkdir test
    cd test/
    wget http://pfam.xfam.org/family/PF00046/hmm    # Downloads the Pfam Homeobox profile HMM
    mv hmm Homeobox.hmm
    python3 ../code/run_hmm_local.py Homeobox.hmm ../database/sequence/pdb_seqres_prot.txt ../database/PDB_mmCIFs/schemes/

The directory `test/HMM_NAME/` contains the color mask JSON files corresponding to the matched domains. These can be loaded into the [LiteMol](https://webchemdev.ncbr.muni.cz/LiteMol/) plugin included in the `html_demo/` directory. The files `test/HMM_NAME.dom` and `test/HMM_NAME.icp` contain the summary of the matched domains and the profile HMM information content profile, respectively. The `HMM_NAME` is read from the contents of the input profile HMM file.

### FASTA input
