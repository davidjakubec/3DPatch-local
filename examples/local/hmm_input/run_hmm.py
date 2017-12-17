from sys import argv
import requests
import subprocess
from os import remove as remove_file

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
                domains.append((chain_ID, hmm_start, hmm_end, seq_start, seq_end, alignments[chain_ID + '/' + seq_start + '-' + seq_end]))
        remove_file('d.out')
        return domains
    domains = process_summary(alignments)
    return domains

def write_domains(domains, hmm_name):
    with open(hmm_name + '.dom', encoding = 'UTF-8', mode = 'w') as f:
        for domain in domains:
            f.write('\t'.join(domain) + '\n')

################################################################################

def main(hmm_file, seq_database_file):
    hmm, hmm_name = read_hmm(hmm_file)
    hmm_information_content_profile = calculate_hmm_information_content_profile(hmm)
    write_hmm_information_content_profile(hmm_information_content_profile, hmm_name)
    domains = run_hmmsearch(hmm_file, seq_database_file)
    write_domains(domains, hmm_name)

if (__name__ == '__main__'):
    main(argv[1], argv[2])
