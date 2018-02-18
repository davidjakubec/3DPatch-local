# 3DPatch-tools

## About

TBD

## Requirements

  - Python >=3.5 (for `subprocess.run()`)
  - [Bio::HMM::Logo](http://skylign.org/help/install) Perl library for calculating HMM information content profiles locally or [Requests](http://docs.python-requests.org/en/master/#) Python 3 library for making [skylign.org](http://skylign.org/) requests
  - [HMMER v3.1b2](http://hmmer.org/) program binaries in `$PATH`
  - `phmmer`/`hmmsearch` target sequence database files; examples can be downloaded using `database/sequence/update_databases.sh`
  - Up-to-date PDB structure residue scheme files; can be generated using `database/PDB_mmCIFs/update_mmCIFs.sh`
  - `run_fasta_local.py` requires `run_hmm_local.py` be present in the same directory; both programs require the `hmm_to_logo.pl` script
  - Write permission in the working directory

## Using 3DPatch-tools

### Profile HMM input

Use the file `code/run_hmm_local.py` as

    python3 run_hmm_local.py PROFILE_HMM_FILE HMMSEARCH_TARGET_DATABASE_FILE STRUCTURE_RESIDUE_SCHEMES_DIRECTORY

The `HMMSEARCH_TARGET_DATABASE_FILE` should be a file containing sequences of PDB chains. Assuming this file and the PDB structure residue schemes have been generated as described in the [Requirements](#requirements), the following example should work

    mkdir test
    cd test/
    wget http://pfam.xfam.org/family/PF00046/hmm    # Downloads the Pfam Homeobox profile HMM
    mv hmm Homeobox.hmm # For clarity, has no effect on naming the output files
    python3 ../code/run_hmm_local.py Homeobox.hmm ../database/sequence/pdb_seqres_prot.txt ../database/PDB_mmCIFs/schemes/

The directory `test/HMM_NAME/` presents for each matched domain a JSON file containing the color mask. These JSON files can be loaded into the [LiteMol](https://webchemdev.ncbr.muni.cz/LiteMol/) plugin included in the `html_demo/` directory. Input profile HMM gathering threshold is used to determine which domain matches are reported. The numbers in the names of the JSON files indicate the first and last residues of the PDB chain sequence included in the alignment to the profile HMM. 

The files `test/HMM_NAME.dom` and `test/HMM_NAME.icp` contain a summary of the matched domains and the profile HMM information content profile, respectively.

`HMM_NAME` (*e.g.*, Homeobox) is read from the contents of the input profile HMM file.

### FASTA input

Use the file `code/run_fasta_local.py` as

    python3 run_fasta_local.py FASTA_FILE PHMMER_TARGET_DATABASE_FILE (HMMSEARCH_TARGET_DATABASE_FILE) STRUCTURE_RESIDUE_SCHEMES_DIRECTORY

`FASTA_FILE` must contain a header and a single protein sequence. The header must be in the format used by the [Protein Data Bank in Europe](https://www.ebi.ac.uk/pdbe/): `>pdb|PDB_ID|CHAIN_ID`. `PDB_ID` and `CHAIN_ID` are used for naming output files.

`PHMMER_TARGET_DATABASE_FILE` should be a file containing high-quality sequences covering the protein sequence space. The script `database/sequence/update_databases.sh` can be used to download the latest release of the UniProtKB/Swiss-Prot database; the UniProt reference proteomes are a good alternative.

`HMMSEARCH_TARGET_DATABASE_FILE` is an optional argument. If it is present, it should be a file containing sequences of PDB chains; `hmmsearch` will be run against these sequences to create mark-up for all PDB structures similar to the FASTA query. If the argument is absent, then the input sequence will be aligned to the profile HMM generated from the `phmmer` search results using `hmmalign`; the mark-up will be generated only for this structure.

If all necessary files have been prepared, the following example should work

    mkdir test2
    cd test2/
    wget http://www.ebi.ac.uk/pdbe/entry/pdb/1ubq/fasta # Downloads FASTA file containing the sequence of ubiquitin
    mv fasta 1ubq_A.fasta # For clarity, has no effect on naming the output files
    python3 ../code/run_fasta_local.py 1ubq_A.fasta ../database/sequence/uniprot_sprot.fasta ../database/sequence/pdb_seqres_prot.txt ../database/PDB_mmCIFs/schemes/

In this case, the directory `test2/PDB_CHAIN_ID/` contains the color mask JSON files corresponding to domains identified using `hmmsearch`. These are the significant hits scoring above per-target and per-domain inclusion thresholds (0.01 and 0.03, respectively).
