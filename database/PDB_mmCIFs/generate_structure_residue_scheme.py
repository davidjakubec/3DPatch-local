from sys import argv
import gzip

################################################################################

def read_CIF_file(CIF_file):
    with gzip.open(CIF_file, mode = 'rt') as f:
        CIF = f.read()
    return CIF

def parse_CIF(CIF):
    blocks = [i.split('\n') for i in CIF.split('\n# \n')]
    structure_residue_scheme = []
    try:
        pdbx_poly_seq_scheme = [i for i in blocks if ((i[0] == 'loop_') and (i[1][:21] == '_pdbx_poly_seq_scheme'))][0]
        pdbx_poly_seq_scheme_headers = [i.strip().split('.')[1] for i in pdbx_poly_seq_scheme if (i[:21] == '_pdbx_poly_seq_scheme')]
        pdbx_poly_seq_scheme_headers_pdb_strand_id_index = pdbx_poly_seq_scheme_headers.index('pdb_strand_id')
        pdbx_poly_seq_scheme_headers_pdb_mon_id_index = pdbx_poly_seq_scheme_headers.index('pdb_mon_id')
        pdbx_poly_seq_scheme_data = [i.strip().split() for i in pdbx_poly_seq_scheme if ((i != 'loop_') and (i[:21] != '_pdbx_poly_seq_scheme'))]
        for line in pdbx_poly_seq_scheme_data:
            structure_residue_scheme.append((line[pdbx_poly_seq_scheme_headers_pdb_strand_id_index], line[pdbx_poly_seq_scheme_headers_pdb_mon_id_index]))
    except IndexError:
        pass
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

def write_structure_residue_scheme(CIF_file, structure_residue_scheme):
    with open('./schemes/' + CIF_file.split('/')[-1].split('.')[0] + '.sch', encoding = 'UTF-8', mode = 'w') as f:
        for residue in structure_residue_scheme:
            f.write('\t'.join(residue) + '\n')

################################################################################

def main(CIF_file):
    print(CIF_file)
    CIF = read_CIF_file(CIF_file)
    structure_residue_scheme = parse_CIF(CIF)
    write_structure_residue_scheme(CIF_file, structure_residue_scheme)

if (__name__ == '__main__'):
    main(argv[1])
