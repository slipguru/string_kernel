import os
import sys
import csv
import pandas as pd

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def shorten(x):
    if len(x) % 3 != 0:
        raise ValueError('Input length should be a multiple of three')

    y = ''
    for i in range(len(x) / 3):
            y += d[x[3 * i:3 * i + 3]]
    return y


def read_pdb(pdb_file, dialect='excel-tab'):
    """Read a protein database file.

    From the pdb file, the heavy and light chain are extracted, with the
    following criteria to divide them in CDRs:
    # LCDR1: 24_L:34_L,
    # LCDR2: 48_L:54_L,
    # LCDR3: 89_L:98_L,
    # HCDR1: 24_H:34_H,
    # HCDR2: 51_H:57_H,
    # HCDR3: 93_H:104_H
    Hence L chain from 24 to 34, from 48 to 54 and so on.
    The numbering is of pdb file, that is the kabat-chothia numbering.
    There can be some IG with ids like 100A, 100B; they must be considered
    if they are inside the intervals specified before.

    Parameters
    ----------
    pdb_file : str
        A protein database file. Delimited according to `dialect`.
    dialect : ('excel-tab', 'excel')
        Dialect for csv.DictReader.

    Returns
    -------
    dict : dictionary
        For each CDR, return the correspondent sequence for the pdb_file
        specified.
    """
    try:
        heavy = []
        light = []
        loops = ['LCDR1', 'LCDR2', 'LCDR3', 'HCDR1', 'HCDR2', 'HCDR3']
        seqs = ['', '', '', '', '', '']
        with open(pdb_file, 'rb') as f:
            for line in f:
                words = line.split("\t")
                if words[0] != 'ATOM':
                    continue
                if words[4] == 'H':
                    if heavy == [] or words[5] != heavy[-1][1]:
                        heavy.append((words[3], words[5]))
                elif words[4] == 'L':
                    if light == [] or words[5] != light[-1][1]:
                        light.append((words[3], words[5]))
                else:
                    raise ValueError
        for pos, t in enumerate(light):
            # remove the suffix A, B, D ...
            n = int(t[1][:-1]) if t[1][-1].isalpha() else int(t[1])
            if 24 <= n <= 34:
                seqs[0] += t[0]
            elif 48 <= n <= 54:
                seqs[1] += t[0]
            elif 89 <= n <= 98:
                seqs[2] += t[0]
        for pos, t in enumerate(heavy):
            # remove the suffix A, B, D ...
            n = int(t[1][:-1]) if t[1][-1].isalpha() else int(t[1])
            if 24 <= n <= 34:
                seqs[3] += t[0]
            elif 51 <= n <= 57:
                seqs[4] += t[0]
            elif 93 <= n <= 104:
                seqs[5] += t[0]

    except IOError:
        sys.exit('ERROR: File %s cannot be read' % pdb_file)
    except Exception as e:
        sys.exit('ERROR: {}'.format(e))
    return {k: v for k, v in zip(loops, map(shorten, seqs))}#, light, heavy

def pdb_to_df(path):

    loops = ['LCDR1', 'LCDR2', 'LCDR3', 'HCDR1', 'HCDR2', 'HCDR3']
    df = pd.DataFrame(columns=loops)

    filenames = [os.path.join(path, f) for f in os.listdir(path)
                 if os.path.isfile(os.path.join(path, f)) and
                 f.endswith('.pdb')]

    for f in filenames:
        dict_seqs = read_pdb(f)
        df = df.append(dict_seqs, ignore_index=True)

    # # To insert the names of the sequences, load them from the merged file
    # # after conversion
    # __df = pd.read_csv('/home/fede/Dropbox/projects/Franco_Fabio_Marcat/'
    #                    'conversioni_ID/tab_final_merged_newid_mutation.csv')
    # indexes = []
    # for f in filenames:
    #     df1[df1['ID TM matrice'] == x]['new_id_tm']
    df.index = map(lambda x: x.split('/')[-1], filenames)
    return df
