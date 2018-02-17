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

#def read_hmm():
#    with open('p.hmm', encoding = 'UTF-8') as f:
#        hmm = f.read()
#    return hmm

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

def run_hmmalign(fasta_file, PDB_chain_ID):
    subprocess.run(['hmmalign', '-o', 'a.sto', '--amino', 'p.hmm', fasta_file])
    remove_file('p.hmm')
    def read_alignment():
        with open('a.sto', encoding = 'UTF-8') as f:
            cont = f.read().split('\n')
        alignment = ''.join([i.split()[-1] for i in cont if (i and (i[:2] != '//') and (i[0] != '#'))])
        model = ''.join([i.split()[-1] for i in cont if (i and (i[:7] == '#=GC RF'))])
        remove_file('a.sto')
        return alignment, model
    alignment, model = read_alignment()
    def process_alignment(alignment, model):
        alignment_position = 0
        model_position = 0
        index_maps = []
        for i in range(len(alignment)):
            if alignment[i] != '-':
                alignment_position += 1
            if model[i] == 'x':
                model_position += 1
                index_maps.append((alignment_position, model_position, i))
        sequence = alignment[index_maps[0][2]:(index_maps[-1][2] + 1)]
        match_state_index_maps = [i for i in index_maps if (alignment[i[2]] != '-')]
        hmm_start, hmm_end, seq_start, seq_end = match_state_index_maps[0][1], match_state_index_maps[-1][1], match_state_index_maps[0][0], match_state_index_maps[-1][0]
        return hmm_start, hmm_end, seq_start, seq_end, sequence
    hmm_start, hmm_end, seq_start, seq_end, sequence = process_alignment(alignment, model)
    domains = [[PDB_chain_ID, str(hmm_start), str(hmm_end), str(seq_start), str(seq_end), sequence]]
    return domains

################################################################################

def main(fasta_file, phmmer_database_file, hmmsearch_database_file, structure_residue_schemes_directory):
    PDB_chain_ID = read_fasta(fasta_file)
    run_phmmer(fasta_file, phmmer_database_file)
    run_hmmbuild()
#    hmm = read_hmm()
#    hmm_information_content_profile = calculate_hmm_information_content_profile(hmm)
    hmm_information_content_profile = calculate_hmm_information_content_profile('p.hmm')
    domains = run_hmmsearch(hmmsearch_database_file)
    calculate_domain_information_content_profiles(domains, hmm_information_content_profile)
    write_domain_color_masks(domains, PDB_chain_ID, structure_residue_schemes_directory)

def alternative_main(fasta_file, phmmer_database_file, structure_residue_schemes_directory):
    PDB_chain_ID = read_fasta(fasta_file)
    run_phmmer(fasta_file, phmmer_database_file)
    run_hmmbuild()
    hmm_information_content_profile = calculate_hmm_information_content_profile('p.hmm')
    domains = run_hmmalign(fasta_file, PDB_chain_ID)
    calculate_domain_information_content_profiles(domains, hmm_information_content_profile)
    write_domain_color_masks(domains, PDB_chain_ID, structure_residue_schemes_directory)

if (__name__ == '__main__'):
    if (len(argv) == 5):
        main(argv[1], argv[2], argv[3], argv[4])
    elif (len(argv) == 4):
        alternative_main(argv[1], argv[2], argv[3])
