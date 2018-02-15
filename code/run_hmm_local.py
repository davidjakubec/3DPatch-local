from sys import argv
from sys import path as sys_path
#import requests
import subprocess
from os import remove as remove_file
from os import mkdir
import json

################################################################################

def read_hmm(hmm_file):
    with open(hmm_file, encoding = 'UTF-8') as f:
        hmm = f.read()
    hmm_name = [i for i in hmm.split('\n') if (i and (i[:4] == 'NAME'))][0].split()[1]
    return hmm, hmm_name

#def calculate_hmm_information_content_profile(hmm):
#    def get_logo_url(hmm):
#        url = 'http://skylign.org/'
#        payload = {'processing': 'hmm'}
#        headers = {'Accept': 'application/json'}
#        files = {'file': ('bla', hmm)}
#        r = requests.post(url, data = payload, headers = headers, files = files)
#        return r.json()['url']
#    url = get_logo_url(hmm)
#    def get_logo(url):
#        headers = {'Accept': 'application/json'}
#        r = requests.get(url, headers = headers)
#        logo = r.json()
#        return logo
#    logo = get_logo(url)
#    hmm_information_content_profile = [round(sum([float(j.split(':')[1]) for j in i]), 3) for i in logo['height_arr']]
#    return hmm_information_content_profile

def calculate_hmm_information_content_profile(hmm_file):
    subprocess.run(['perl', sys_path[0] + '/hmm_to_logo.pl', hmm_file])
    def read_logo():
        with open('logo.json', encoding = 'UTF-8') as f:
            logo = json.load(f)
        remove_file('logo.json')
        return logo
    logo = read_logo()
    hmm_information_content_profile = [round(sum([float(j.split(':')[1]) for j in i]), 3) for i in logo['height_arr']]
    return hmm_information_content_profile

def write_hmm_information_content_profile(hmm_information_content_profile, hmm_name):
    with open(hmm_name + '.icp', encoding = 'UTF-8', mode = 'w') as f:
        for i in hmm_information_content_profile:
            f.write(str(i) + '\n')

def run_hmmsearch(hmm_file, seq_database_file):
    subprocess.run(['hmmsearch', '-o', '/dev/null', '-A', 'a.sto', '--domtblout', 'd.out', '--notextw', '--cut_ga', hmm_file, seq_database_file])
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

def write_domains(domains, hmm_name):
    with open(hmm_name + '.dom', encoding = 'UTF-8', mode = 'w') as f:
        for domain in domains:
            f.write('\t'.join(domain) + '\n')

def calculate_domain_information_content_profiles(domains, hmm_information_content_profile):
    for domain in domains:
        seq = domain[5]
        domain_information_content_profile = []
        match_state_count = 0
        for letter in seq:
            if (letter == '-'):
                match_state_count += 1
            elif (letter.isupper() == True):
                domain_information_content_profile.append(hmm_information_content_profile[match_state_count])
                match_state_count += 1
            elif (letter.islower() == True):
                domain_information_content_profile.append('i')
        domain.append(domain_information_content_profile)

def write_domain_color_masks(domains, hmm_name, structure_residue_schemes_directory):
    try:
        mkdir(hmm_name)
    except FileExistsError:
        pass
    def read_structure_residue_scheme(PDB_ID, structure_residue_schemes_directory):
        with open(structure_residue_schemes_directory + PDB_ID + '.sch', encoding = 'UTF-8') as f:
            structure_residue_scheme = [i.split('\t') for i in f.read().split('\n') if i]
        return structure_residue_scheme
    def calculate_structure_residue_scheme_information_content_profile(chain_ID, seq_start, seq_end, domain_information_content_profile, structure_residue_scheme):
        structure_residue_scheme_information_content_profile = []
        chain_residue_index = 1
        domain_residue_index = 0
        for residue in structure_residue_scheme:
            if (residue[0] != chain_ID):
                structure_residue_scheme_information_content_profile.append('d')
            else:
                if ((chain_residue_index >= seq_start) and (chain_residue_index <= seq_end)):
                    structure_residue_scheme_information_content_profile.append(domain_information_content_profile[domain_residue_index])
                    domain_residue_index += 1
                else:
                    structure_residue_scheme_information_content_profile.append('m')
                chain_residue_index += 1
        return structure_residue_scheme_information_content_profile
    def convert_profile_to_color_mask(structure_information_content_profile):
        domain_color_mask = ['skipped',]
        for information_content in structure_information_content_profile:
            if (information_content == 'd'):
                color = 'd'
            elif (information_content == 'm'):
                color = 'm'
            elif (information_content == 'i'):
                color = 'i'
            else:
                normalized_information_content = information_content / 6.45311498641968
                color = str(int(10 * normalized_information_content))
            domain_color_mask.append(color)
        return domain_color_mask
    for domain in domains:
        PDB_ID = domain[0].split('_')[0]
        chain_ID = domain[0].split('_')[1]
        seq_start = int(domain[3])
        seq_end = int(domain[4])
        domain_information_content_profile = domain[6]
        structure_residue_scheme = read_structure_residue_scheme(PDB_ID, structure_residue_schemes_directory)
        structure_residue_scheme_information_content_profile = calculate_structure_residue_scheme_information_content_profile(chain_ID, seq_start, seq_end, domain_information_content_profile, structure_residue_scheme)
        structure_information_content_profile = [structure_residue_scheme_information_content_profile[i] for i in range(len(structure_residue_scheme_information_content_profile)) if (structure_residue_scheme[i][1] != '?')]
        domain_color_mask = convert_profile_to_color_mask(structure_information_content_profile)
        with open(hmm_name + '/' + PDB_ID + '_' + chain_ID + '_' + str(seq_start) + '-' + str(seq_end) + '.json', encoding = 'UTF-8', mode = 'w') as f:
            json.dump({'pdbId': PDB_ID, 'colorMask': domain_color_mask}, f)

################################################################################

def main(hmm_file, seq_database_file, structure_residue_schemes_directory):
    hmm, hmm_name = read_hmm(hmm_file)
#    hmm_information_content_profile = calculate_hmm_information_content_profile(hmm)
    hmm_information_content_profile = calculate_hmm_information_content_profile(hmm_file)
    write_hmm_information_content_profile(hmm_information_content_profile, hmm_name)
    domains = run_hmmsearch(hmm_file, seq_database_file)
    write_domains(domains, hmm_name)
    calculate_domain_information_content_profiles(domains, hmm_information_content_profile)
    write_domain_color_masks(domains, hmm_name, structure_residue_schemes_directory)

if (__name__ == '__main__'):
    main(argv[1], argv[2], argv[3])
