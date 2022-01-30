import re
import numpy as np
import matplotlib.pyplot as plt

dict_aminoacid_composition = {  # format = (C,H,N,O,S,C13,D) and be edited to create modified amino acids
    "A": (3, 5, 1, 1, 0),
    "R": (6, 12, 4, 1, 0),
    "N": (4, 6, 2, 2, 0),
    "D": (4, 5, 1, 3, 0),
    "C": (3, 5, 1, 1, 1),
    "E": (5, 7, 1, 3, 0),
    "Q": (5, 8, 2, 2, 0),
    "G": (2, 3, 1, 1, 0),
    "H": (6, 7, 3, 1, 0),
    "I": (6, 11, 1, 1, 0),
    "L": (6, 11, 1, 1, 0),
    "K": (6, 12, 2, 1, 0),
    "M": (5, 9, 1, 1, 1),
    "F": (9, 9, 1, 1, 0),
    "P": (5, 7, 1, 1, 0),
    "S": (3, 5, 1, 2, 0),
    "T": (4, 7, 1, 2, 0),
    "U": (5, 5, 1, 2, 0),
    "W": (11, 10, 2, 1,0),
    "Y": (9, 9, 1, 2, 0),
    "V": (5, 9, 1, 1, 0)
}
dict_aminoacid_aveMasses = {  # format = (C,H,N,O,S,C13,D) and be edited to create modified amino acids
    "A": 71.0779,
    "R": 156.1857,
    "N": 114.1026,
    "D": 115.0874,
    "C": 103.1429,
    "E": 129.114,
    "Q": 128.1292	,
    "G": 57.0513,
    "H": 137.1393,
    "I": 113.1576	,
    "L": 113.1576	,
    "K": 128.1723,
    "M": 131.1961,
    "F": 147.1739,
    "P": 97.1152,
    "S": 87.0773,
    "T": 	101.1039,
    "U": 150.0379	,
    "W": 186.2099,
    "Y": 163.1733,
    "V": 	99.1311}
dict_atomic_mass = {
    # mono isotopic mass for each atom - can expand for isotopes
    "C": 12,
    "H": 1.007825,
    "N": 14.003074,
    "O": 15.994915,
    "S": 31.972072}
expasy_rules = {
    'arg-c':         r'R',
    'asp-n':         r'\w(?=D)',
    'bnps-skatole' : r'W',
    'caspase 1':     r'(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
    'caspase 2':     r'(?<=DVA)D(?=[^PEDQKR])',
    'caspase 3':     r'(?<=DMQ)D(?=[^PEDQKR])',
    'caspase 4':     r'(?<=LEV)D(?=[^PEDQKR])',
    'caspase 5':     r'(?<=[LW]EH)D',
    'caspase 6':     r'(?<=VE[HI])D(?=[^PEDQKR])',
    'caspase 7':     r'(?<=DEV)D(?=[^PEDQKR])',
    'caspase 8':     r'(?<=[IL]ET)D(?=[^PEDQKR])',
    'caspase 9':     r'(?<=LEH)D',
    'caspase 10':    r'(?<=IEA)D',
    'chymotrypsin high specificity' : r'([FY](?=[^P]))|(W(?=[^MP]))',
    'chymotrypsin low specificity':
        r'([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
    'clostripain':   r'R',
    'cnbr':          r'M',
    'enterokinase':  r'(?<=[DE]{3})K',
    'factor xa':     r'(?<=[AFGILTVM][DE]G)R',
    'formic acid':   r'D',
    'glutamyl endopeptidase': r'E',
    'granzyme b':    r'(?<=IEP)D',
    'hydroxylamine': r'N(?=G)',
    'iodosobenzoic acid': r'W',
    'lysc':          r'K',
    'ntcb':          r'\w(?=C)',
    'pepsin ph1.3':  r'((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
                     r'((?<=[^HKR][^P])[FL](?=\w[^P]))',
    'pepsin ph2.0':  r'((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
                     r'((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
    'proline endopeptidase': r'(?<=[HKR])P(?=[^P])',
    'proteinase k':  r'[AEFILTVWY]',
    'staphylococcal peptidase i': r'(?<=[^E])E',
    'thermolysin':   r'[^DE](?=[AFILMV])',
    'thrombin':      r'((?<=G)R(?=G))|'
                     r'((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin':       r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
    'trypsin_exception': r'((?<=[CD])K(?=D))|((?<=C)K(?=[HY]))|((?<=C)R(?=K))|((?<=R)R(?=[HR]))',
    }




def Create_Fasta_dict(file_name):
    dict_Fasta = {}
    f = open(file_name,"r")
    Fasta_file = f.read().split("\n>")
    for chunk in Fasta_file:
        lines=chunk.split('\n')
        if lines[0].startswith('Rev'):
            ""
        else:
            try:
                uniprot_name=lines[0].split("|")[1]
            except:
                uniprot_name=lines[0]
            seq="".join(lines[1:])

            dict_Fasta[uniprot_name]=seq
    f.close()

    return dict_Fasta

def generate_pattern_dict(max_length = 50):
    dict_forward_coordinates = {}
    for j in range(1, max_length):
        list_len3 = []
        for i in range(1, j + 1):
            list_len3.append([0, i])
        dict_forward_coordinates[j] = list_len3
    return dict_forward_coordinates


def peptide_generator(seq,enzyme = "trypsin",missed = 0,min_len = 0,max_len = 25,max_length = 50):
    cut_re_expasy = expasy_rules[enzyme]
    dict_forward_coordinates = generate_pattern_dict(max_length)

    peptide_list = [m.end() for m in re.finditer(cut_re_expasy, seq)]
    peptide_list.insert(0, 0)
    peptide_list.append(len(seq))

    list_coordinates_forward = [(x[0] + peptide_list[i], x[1] + peptide_list[i]) for j in range(1, missed + 2) for i in
                                range(1, len(peptide_list) - j) if peptide_list[i + j] - peptide_list[i] < max_length
                                for x in dict_forward_coordinates[peptide_list[i + j] - peptide_list[i]] if
                                min_len < (x[1] - x[0]) < max_len]
    list_coordinates_backwards = [(peptide_list[i] - x[1], peptide_list[i] - x[0]) for j in range(1, missed + 2) for i
                                  in range(1, len(peptide_list) - 1) if i - j >= 0 if
                                  peptide_list[i] - peptide_list[i - j] < max_length for x in
                                  dict_forward_coordinates[peptide_list[i] - peptide_list[i - j]] if
                                  peptide_list[i] - peptide_list[i - j] <= max_length if
                                  min_len < (x[1] - x[0]) < max_len]

    list_total_coordinates = list_coordinates_forward + list_coordinates_backwards

    list_peptides = [seq[list_total_coordinates[i][0]:list_total_coordinates[i][1]] for i in
                     range(len(list_total_coordinates))]

    return set(list_peptides)



def create_reverse_fasta(file_name):   #generates a dictonary from a fasta file
    dict_fasta = {}
    f = open(file_name,"r")     #opens the file

    fasta_file = f.read().split("\n>")      #reads the file and splits into chunks - \n> is before every entry

    for chunk in fasta_file:        # go though chunk by chunk
        lines=chunk.split('\n')         #split by lines - first line contains protein name and other junk
        try:
            uniprot_name=lines[0].split(">")[1]     # catches the first one which will have extra ">" in it otherwise
        except:
            uniprot_name=lines[0]       # if there is no | - just use the whole line
        seq = "".join(lines[1:])      # join all but the first line together to give the sequence
        dict_fasta[uniprot_name]= seq       # create dictionary
    f.close()
    write_file = open("rev_" +file_name,"w")
    for each in dict_fasta:
        write_file.write(">" + each + "\n")
        write_file.write(dict_fasta[each] + "\n")
        write_file.write(">rev_" + each + "\n")
        write_file.write(dict_fasta[each][::-1] + "\n")
    write_file.close
    return

def Create_dict_pepXML(file_name):
    dict_xml_peptides = {}
    xml_file = open(file_name,"r")
    data = xml_file.read().split('peptide="')
    for i in range(1,len(data)):
        peptide_seq = data[i].split('"')[0]
        protein_ID =  data[i].split('"')[20].split("|")[1]
        mass_dif = data[i].split('"')[2]
        calc_neutral_pep_mass = data[i].split('"')[4]
        missed = data[i].split('"')[8]
        num_matched_ions = data[i].split('"')[18]

        dict_xml_peptides[peptide_seq] = {}
        dict_xml_peptides[peptide_seq]["Protein_ID"]= protein_ID

        dict_xml_peptides[peptide_seq]["Mass_diff"] = mass_dif
        dict_xml_peptides[peptide_seq]["Calc_pep_mass"] = calc_neutral_pep_mass
        dict_xml_peptides[peptide_seq]["missed"] = missed
        dict_xml_peptides[peptide_seq]["number_matched_ions"] = num_matched_ions
        if protein_ID.split("|")[0][0:3] != "rev":
            dict_xml_peptides[peptide_seq]["reverse_seq_match"] = False
        else:
            dict_xml_peptides[peptide_seq]["reverse_seq_match"] = True
    return dict_xml_peptides


def peak_picker(X, Y, view_plot=0):
    xvals = np.linspace(X[0], X[-1], 1000)
    Y_interp = np.interp(xvals, X, Y)
    Y_1stD = np.gradient(Y_interp)
    peaks = []
    for i in range(len(Y_1stD) - 1):
        if Y_1stD[i] >= 0.1 and Y_1stD[i + 1] <= 0.1:
            peaks.append((xvals[i], Y_interp[i]))

    if view_plot == 1:
        plt.plot(X, Y, label="data")

        for i in range(len(peaks)):
            plt.plot([peaks[i][0]] * 2, [peaks[i][1], peaks[i][1]], 'b-x')
        plt.legend()
        plt.show()

    return peaks




def mono_mass_calc(seq, charge=1):
    """
    Returns the monoisotopic mass for a peptide sequence in the form GAP(15.994915)GP(15.994915)DGNNG

    mass modifications can be placed after the amino acid code that is modified inside brackets (XXX.XXXXXX) up to 6 d.p.
    and to a max mass of 999.999999.For positive modifications a + symbol is not required
    however for negative modifications a - symbol must be before the mass

    The charge can be specified or will default as +1

    if the input sequence contains the previous and next residues (such as R.GAPQGVQGGK.G) they will be ignored

"""

    try:
        if seq[1] == ".":
            seq = seq[2:-2]
    except:
        pass
    peptide_comp = np.array([0, 2, 0, 1, 0])
    for i in dict_aminoacid_composition:
        for match in re.finditer(i, seq):
            peptide_comp += dict_aminoacid_composition[i]
    v = peptide_comp[0]  # C
    w = peptide_comp[1]  # H
    x = peptide_comp[2]  # N
    y = peptide_comp[3]  # O
    z = peptide_comp[4]  # S

    main_mass_peak = v * dict_atomic_mass["C"] + w * dict_atomic_mass["H"] + x * \
                     dict_atomic_mass["N"] + y * dict_atomic_mass["O"] + z * \
                     dict_atomic_mass["S"]

    for match in re.finditer('[+-]?(\d{1,3}\.\d{1,6})', seq):
        main_mass_peak += float(match.group())
    main_mass_peak = float(((main_mass_peak + charge * 1.007276466) / charge))

    return main_mass_peak


def Y_B_ion_generator(seq,max_charge=3):

    if seq[1] == ".":
        seq = seq[2:-2]
    b_peptides = [seq[:i+1]for i in range(0,len(seq)-1) if seq[i+1] in dict_aminoacid_composition]
    b_peptides.append(seq)

    dict_ions = {}
    dict_ions["B"] = {}
    for i in range(1,max_charge):
        dict_ions["B"][i] = {}
        for j in range(0,len(b_peptides)):
            dict_ions["B"][i][j+1] = mono_mass_calc(b_peptides[j],i) - (15.994915 + 2*1.007825)

    dict_ions["Y"] = {}
    y_peptides = [seq[i-1:]for i in range(1,len(seq)+1) if seq[i-1] in dict_aminoacid_composition ]
    for i in range(1,max_charge):
        dict_ions["Y"][i] = {}
        for j in range(0,len(y_peptides)):
            dict_ions["Y"][i][len(y_peptides)-j] = mono_mass_calc(y_peptides[j],i)
    return dict_ions