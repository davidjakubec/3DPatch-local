from sys import argv
import requests
import subprocess
from os import remove as remove_file
from os import mkdir
import json

################################################################################

def read_hmm(hmm_file):
    with open(hmm_file, encoding = 'UTF-8') as f:
        hmm = f.read();
    hmm_name = [i for i in hmm.split('\n') if (i and (i[:4] == 'NAME'))][0].split()[1]
    return hmm, hmm_name

def calculate_hmm_information_content_profile(hmm):
    def get_logo_url(hmm):
        url = 'http://skylign.org/'
        payload = {'processing': 'hmm'}
        headers = {'Accept': 'application/json'}
        files = {'file': ('bla', hmm)}
        r = requests.post(url, data = payload, headers = headers, files = files)
        return r.json()['url']
    url = get_logo_url(hmm)
    def get_logo(url):
        headers = {'Accept': 'application/json'}
        r = requests.get(url, headers = headers)
        logo = r.json()
        return logo
    logo = get_logo(url)
    hmm_information_content_profile = [round(sum([float(j.split(':')[1]) for j in i]), 3) for i in logo['height_arr']]
    return hmm_information_content_profile

def write_hmm_information_content_profile(hmm_information_content_profile, hmm_name):
    with open(hmm_name + '.icp', encoding = 'UTF-8', mode = 'w') as f:
        f.write(' '.join([str(i) for i in hmm_information_content_profile]))

def run_hmmsearch(hmm_file, seq_database_file):
    subprocess.run(['hmmsearch', '-o', '/dev/null', '-A', 'a.sto', '--domtblout', 'd.out', hmm_file, seq_database_file])
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

def read_color_scale(color_scale_file):
    with open(color_scale_file, encoding = 'UTF-8') as f:
        color_scale = [[float(j) for j in i.split()] for i in f.read().split('\n') if i]
    return color_scale

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

def download_and_parse_CIFs(domains):
    PDB_IDs = set([i[0].split('_')[0] for i in domains])
    def download_CIF(PDB_ID):
        url = 'https://www.ebi.ac.uk/pdbe/static/entry/' + PDB_ID + '_updated.cif'
        r = requests.get(url)
        return r.text
    def parse_CIF(CIF):
        blocks = [i.split('\n') for i in CIF.split('\n#\n')]
        pdbx_poly_seq_scheme = [i for i in blocks if ((i[0] == 'loop_') and (i[1][:21] == '_pdbx_poly_seq_scheme'))][0]
        pdbx_poly_seq_scheme_headers = [i.strip().split('.')[1] for i in pdbx_poly_seq_scheme if (i[:21] == '_pdbx_poly_seq_scheme')]
        pdbx_poly_seq_scheme_headers_pdb_strand_id_index = pdbx_poly_seq_scheme_headers.index('pdb_strand_id')
        pdbx_poly_seq_scheme_headers_pdb_mon_id_index = pdbx_poly_seq_scheme_headers.index('pdb_mon_id')
        pdbx_poly_seq_scheme_data = [i.strip().split() for i in pdbx_poly_seq_scheme if ((i != 'loop_') and (i[:21] != '_pdbx_poly_seq_scheme'))]
        structure_residue_scheme = []
        for line in pdbx_poly_seq_scheme_data:
            structure_residue_scheme.append((line[pdbx_poly_seq_scheme_headers_pdb_strand_id_index], line[pdbx_poly_seq_scheme_headers_pdb_mon_id_index]))
        try:
            pdbx_nonpoly_scheme = [i for i in blocks if ((i[0] == 'loop_') and (i[1][:20] == '_pdbx_nonpoly_scheme'))][0]
            pdbx_nonpoly_scheme_headers = [i.strip().split('.')[1] for i in pdbx_nonpoly_scheme if (i[:20] == '_pdbx_nonpoly_scheme')]
            pdbx_nonpoly_scheme_headers_pdb_strand_id_index = pdbx_nonpoly_scheme_headers.index('pdb_strand_id')
            pdbx_nonpoly_scheme_headers_pdb_mon_id_index = pdbx_nonpoly_scheme_headers.index('pdb_mon_id')
            pdbx_nonpoly_scheme_data = [i.strip().split() for i in pdbx_nonpoly_scheme if ((i != 'loop_') and (i[:20] != '_pdbx_nonpoly_scheme'))]
            for line in pdbx_nonpoly_scheme_data:
                structure_residue_scheme.append((line[pdbx_nonpoly_scheme_headers_pdb_strand_id_index], line[pdbx_nonpoly_scheme_headers_pdb_mon_id_index]))
        except IndexError:
            pass
        return structure_residue_scheme
    structure_residue_schemes = dict()
    for PDB_ID in PDB_IDs:
        CIF = download_CIF(PDB_ID)
        structure_residue_scheme = parse_CIF(CIF)
        structure_residue_schemes[PDB_ID] = structure_residue_scheme
    return structure_residue_schemes

def calculate_domain_color_masks(domains, structure_residue_schemes, color_scale):
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
    def convert_profile_to_colors(structure_information_content_profile, color_scale):
        domain_color_mask = ['skipped',]
        for information_content in structure_information_content_profile:
            if (information_content == 'd'):
                color = {'r': 0.839, 'g': 0.910, 'b': 0.976}
            elif (information_content == 'm'):
                color = {'r': 0.749, 'g': 0.937, 'b': 0.561}
            elif (information_content == 'i'):
                color = {'r': 0.749, 'g': 0.937, 'b': 0.561}
            else:
                normalized_information_content = information_content / 6.45311498641968
                color_index = round(255 * normalized_information_content)
                color = {'r': color_scale[color_index][0], 'g': color_scale[color_index][1], 'b': color_scale[color_index][2]}
            domain_color_mask.append(color)
        return domain_color_mask
    for domain in domains:
        PDB_ID = domain[0].split('_')[0]
        chain_ID = domain[0].split('_')[1]
        seq_start = int(domain[3])
        seq_end = int(domain[4])
        domain_information_content_profile = domain[6]
        structure_residue_scheme = structure_residue_schemes[PDB_ID]
        structure_residue_scheme_information_content_profile = calculate_structure_residue_scheme_information_content_profile(chain_ID, seq_start, seq_end, domain_information_content_profile, structure_residue_scheme)
        structure_information_content_profile = [structure_residue_scheme_information_content_profile[i] for i in range(len(structure_residue_scheme_information_content_profile)) if (structure_residue_scheme[i][1] != '?')]
        domain_color_mask = convert_profile_to_colors(structure_information_content_profile, color_scale)
        domain.append(domain_color_mask)

def write_domain_color_masks(domains, hmm_name):
    try:
        mkdir(hmm_name)
    except FileExistsError:
        pass
    for domain in domains:
        chain_ID = domain[0]
        seq_start = domain[3]
        seq_end = domain[4]
        domain_color_mask = domain[7]
        with open(hmm_name + '/' + chain_ID + '_' + seq_start + '-' + seq_end + '.json', encoding = 'UTF-8', mode = 'w') as f:
            json.dump({'pdbId': chain_ID.split('_')[0], 'colorMask': domain_color_mask}, f)

################################################################################

def main(hmm_file, seq_database_file, color_scale_file):
    hmm, hmm_name = read_hmm(hmm_file)
    hmm_information_content_profile = calculate_hmm_information_content_profile(hmm)
    write_hmm_information_content_profile(hmm_information_content_profile, hmm_name)
    domains = run_hmmsearch(hmm_file, seq_database_file)
    write_domains(domains, hmm_name)
    calculate_domain_information_content_profiles(domains, hmm_information_content_profile)
    structure_residue_schemes = download_and_parse_CIFs(domains)
    color_scale = read_color_scale(color_scale_file)
    calculate_domain_color_masks(domains, structure_residue_schemes, color_scale)
    write_domain_color_masks(domains, hmm_name)

if (__name__ == '__main__'):
    main(argv[1], argv[2], argv[3])
