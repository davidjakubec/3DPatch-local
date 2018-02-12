from sys import argv
import subprocess
from os import remove as remove_file
from run_hmm_local import calculate_hmm_information_content_profile
from run_hmm_local import calculate_domain_information_content_profiles
from run_hmm_local import write_domain_color_masks

################################################################################

def read_fasta(fasta_file):
    with open(fasta_file, encoding = 'UTF-8') as f:
        cont = [i for i in f.read().split('\n') if i]
    PDB_chain_ID = '_'.join(cont[0].split('|')[-2:])
#    chain_sequence = ''.join(cont[1:])
    return PDB_chain_ID

def run_phmmer(fasta_file, phmmer_database_file):
    subprocess.run(['phmmer', '-o', '/dev/null', '-A', 'a.sto', '--notextw', '--incE', '0.01', '--incdomE', '0.03', fasta_file, phmmer_database_file])

def run_hmmbuild():
    subprocess.run(['hmmbuild', '-o', '/dev/null', '--amino', '--wpb', '--eent', 'p.hmm', 'a.sto'])
    remove_file('a.sto')

def read_hmm():
    with open('p.hmm', encoding = 'UTF-8') as f:
        hmm = f.read()
    return hmm

def run_hmmsearch(hmmsearch_database_file):
    subprocess.run(['hmmsearch', '-o', '/dev/null', '-A', 'a.sto', '--domtblout', 'd.out', '--notextw', '--incE', '0.01', '--incdomE', '0.03', 'p.hmm', hmmsearch_database_file])
    remove_file('p.hmm')
    def process_alignments():
        with open('a.sto', encoding = 'UTF-8') as f:
            cont = [i.split() for i in f.read().split('\n') if (i and (i[0] != '#') and (i[:2] != '//'))]
        alignments = dict()
        for line in cont:
            alignments[line[0]] = line[1]
        remove_file('a.sto')
        return alignments
    alignments = process_alignments()
    def process_summary(alignments):
        with open('d.out', encoding = 'UTF-8') as f:
            cont = [i.split() for i in f.read().split('\n') if (i and (i[0] != '#'))]
        domains = []
        for line in cont:
            chain_ID = line[0]
            hmm_start, hmm_end, seq_start, seq_end = line[15:19]
            if ((chain_ID + '/' + seq_start + '-' + seq_end) in alignments.keys()):
                domains.append([chain_ID, hmm_start, hmm_end, seq_start, seq_end, alignments[chain_ID + '/' + seq_start + '-' + seq_end]])
        remove_file('d.out')
        return domains
    domains = process_summary(alignments)
    return domains

################################################################################

def main(fasta_file, phmmer_database_file, hmmsearch_database_file, structure_residue_schemes_directory):
    PDB_chain_ID = read_fasta(fasta_file)
    run_phmmer(fasta_file, phmmer_database_file)
    run_hmmbuild()
    hmm = read_hmm()
    hmm_information_content_profile = calculate_hmm_information_content_profile(hmm)
    domains = run_hmmsearch(hmmsearch_database_file)
    calculate_domain_information_content_profiles(domains, hmm_information_content_profile)
    write_domain_color_masks(domains, PDB_chain_ID, structure_residue_schemes_directory)

if (__name__ == '__main__'):
    main(argv[1], argv[2], argv[3], argv[4])
