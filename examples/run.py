from sys import argv
import logging
import requests

############################################################

def parse_fasta(fasta_file):
    with open(fasta_file, encoding = 'UTF-8') as f:
        cont = [i for i in f.read().split('\n') if i]
    fasta_header = cont[0].lstrip(">")
    seq = ''.join(cont[1:])
    return fasta_header, seq

def configure_logging(fasta_header):
    logging.basicConfig(filename = fasta_header + '.log', filemode = 'w', level = logging.INFO)

def phmmer_search(seq):
    url = 'https://www.ebi.ac.uk/Tools/hmmer/search/phmmer'
    payload = {'algo': 'phmmer', 'seq': seq, 'seqdb': 'uniprotrefprot'}
    headers = {'Accept': 'application/json'}
    r = requests.post(url, data = payload, headers = headers, allow_redirects = False)
    phmmer_search_job_id = r.headers['Location'].split('/')[6]
    return phmmer_search_job_id

def hmmsearch_search(phmmer_search_job_id):
    url = 'https://www.ebi.ac.uk/Tools/hmmer//search/hmmsearch?uuid=' + phmmer_search_job_id + '.1'
    payload = {'algo': 'hmmsearch', 'seqdb': 'pdb'}
    headers = {'Accept': 'application/json'}
    r = requests.post(url, data = payload, headers = headers)
    hmmsearch_search_hits = r.json()['results']['hits']



############################################################

def main(fasta_file):
    fasta_header, seq = parse_fasta(fasta_file)
    configure_logging(fasta_header)
    phmmer_search_job_id = phmmer_search(seq)
    hmmsearch_search(phmmer_search_job_id)


if (__name__ == "__main__"):
    main(argv[1])
