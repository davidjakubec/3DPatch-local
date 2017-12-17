from sys import argv
import requests

################################################################################

def read_hmm(hmm_file):
    with open(hmm_file, encoding = 'UTF-8') as f:
        hmm = f.read()
    return hmm

def calculate_hmm_information_content_profile(hmm):
    def get_logo_url(hmm):
        url = 'http://skylign.org/'
        payload = {'processing': 'hmm'}
        headers = {'Accept': 'application/json'}
        files = {'file': ('run-fasta.hmm', hmm)}
        r = requests.post(url, data = payload, headers = headers, files = files)
        return r.json()['url']
    url = get_logo_url(hmm)
    def get_logo(url):
        headers = {'Accept': 'application/json'}
        r = requests.get(url, headers = headers)
        logo = r.json()
        return logo
    logo = get_logo(url)
    hmm_information_content_profile = [sum([float(j.split(':')[1]) for j in i]) for i in logo['height_arr']]
    return hmm_information_content_profile



################################################################################

def main(hmm_file):
    hmm = read_hmm(hmm_file)
    hmm_information_content_profile = calculate_hmm_information_content_profile(hmm)



if (__name__ == "__main__"):
    main(argv[1])
