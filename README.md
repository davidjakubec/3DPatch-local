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

The directory `test/HMM_NAME/` contains the color mask JSON files corresponding to the matched domains. These domains are the significant hits scoring above the profile HMM gathering threshold. The numbers in the names of these files indicate the first and last residues of the PDB chain sequence included in the alignment to the profile HMM. These files can be loaded into the [LiteMol](https://webchemdev.ncbr.muni.cz/LiteMol/) plugin included in the `html_demo/` directory. The files `test/HMM_NAME.dom` and `test/HMM_NAME.icp` contain the summary of the matched domains and the profile HMM information content profile, respectively. The `HMM_NAME` is read from the contents of the input profile HMM file.

### FASTA input

Use the file `code/run_fasta_local.py` as

    python3 run_fasta_local.py FASTA_FILE PHMMER_TARGET_DATABASE_FILE HMMSEARCH_TARGET_DATABASE_FILE STRUCTURE_RESIDUE_SCHEMES_DIRECTORY

The `PHMMER_TARGET_DATABASE_FILE` should be a file containing high-quality sequences covering the protein sequence space. The script `database/sequence/update_databases.sh` can be used to download the latest release of the UniProtKB/Swiss-Prot database; the UniProt reference proteomes are a good alternative. `FASTA_FILE` must contain a header and a single protein sequence. The header must be in the format used by the [Protein Data Bank in Europe](https://www.ebi.ac.uk/pdbe/): `>pdb|PDB_ID|CHAIN_ID`; `PDB_ID` and `CHAIN_ID` are used for naming the output. Assuming all necessary files are present, the following example should work

    mkdir test2
    cd test2/
    wget http://www.ebi.ac.uk/pdbe/entry/pdb/1ubq/fasta # Downloads FASTA file containing the sequence of ubiquitin
    mv fasta 1ubq_A.fasta
    python3 ../code/run_fasta_local.py 1ubq_A.fasta ../database/sequence/uniprot_sprot.fasta ../database/sequence/pdb_seqres_prot.txt ../database/PDB_mmCIFs/schemes/

The directory `test2/PDB_CHAIN_ID/` contains the color mask JSON files corresponding to the domains identified using `hmmsearch`. These are the significant hits scoring above the per-target and per-domain inclusion thresholds (0.01 and 0.03, respectively).
