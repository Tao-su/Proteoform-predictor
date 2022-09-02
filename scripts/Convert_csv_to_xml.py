'''
This script is used to convert csv file downloaded from ATCC and containing aa sequence to Uniprot XML file and fasta file for PTM_predictor0
'''
import xml.etree.ElementTree as ET
import hashlib
import proteomic_tools as PT
from xml.dom import minidom
import pandas as pd
import math
file = "ATCC13880_hypothetical_proteins.xml"
ATCC_protein_dataframe = pd.read_csv(r'D:\Kelleher_Lab\005_PTM_predictor\ATCC_13880.csv')
dummy_index = 0
dict = {}
for i in range(ATCC_protein_dataframe.shape[0]):
    if not pd.notna(ATCC_protein_dataframe['Uniprot ID'][i]) and pd.notna(ATCC_protein_dataframe['Amino Acid Sequence'][i]):
        dummy_index += 1
        dummy_acc = 'O' + str(dummy_index).zfill(5)
        dict[ATCC_protein_dataframe['Product'][i]+'_'+dummy_acc] = ATCC_protein_dataframe['Amino Acid Sequence'][i]

def GenerateXML(file='test.xml',
                accession='HC_AAACGGGAGATCCTGT-1_contig_2_IGHD_IGHD*02',
                name='HC_AAACGGGAGATCCTGT-1_contig_2_IGHD_IGHD*02',
                fullName='SequnceID=2 Patient=H1 C=IGHD V=IGHV4-34 D=None J=IGHJ2',
                sequence='MDLLHKNMKHLWFFLLLVAAPRWVLSQVQLQQWGAGLLKPSETLSLTCAVYGGSFSGYYWSWIRQPPGKGLEWIGEINHSGSTNYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCARGGGYGERNWYFDLWGRGTLVTVSSAPTKAPDVFPIISGCRHPKDNSPVVLACLITGYHPTSVTVTWYMGTQSQPQRTFPEIQRRDSYYMTSSQLSTPLQQWRQGEYKCVVQHTASKSKKEIFRWPESPKAQASSVPTAQPQAEGSLAKATTAPATTRNTGRGGEEKKKEKEKEEQEERETKTPECPSHTQPLGVYLLTPAVQDLWLRDKATFTCFVVGSDLKDAHLTWEVAGKVPTGGVEEGLLERHSNGSQSQHSRLTLPRSLWNAGTSVTCTLNHPSLPPQRLMALREPAAQAPVKLSLNLLASSDPPEAASWLLCEVSGFSPPNILLMWLEDQREVNTSGFAPARPPPQPRSTTFWAWSVLRVPAPPSPQPATYTCVVSHEDSRTLLNASRSLEVSYVTDHGPMK'):
    '''
    This function takes protein sequence information as input to generate a Uniprot xml root which is appended to previous xml file.
    Parameters
    ----------
    file: location of the new file or previously generated file
    accession: acc of the protein
    name: name of the protein
    fullName: full name of the protein
    sequence: sequence of the protein

    Returns: no return value but an xml file will be written
    -------

    '''
    root = ET.Element("entry", dataset="Swiss-Prot", created="2022-01-16", modified="2022-01-16", version="1")
    l1_1 = ET.Element("accession")
    root.append(l1_1)
    l1_1.text = accession

    l1_2 = ET.Element("name")
    root.append(l1_2)
    l1_2.text = name

    l1_3 = ET.Element("protein")
    root.append(l1_3)
    l2_1 = ET.SubElement(l1_3, "recommendedName")
    l3_1 = ET.SubElement(l2_1, "fullName")
    l3_1.text = fullName

    l1_4 = ET.Element("sequence", length=str(len(sequence)), mass=str(int(math.ceil(PT.mono_mass_calc(sequence)))), checksum=hashlib.md5(sequence.encode('utf-8')).hexdigest(), modified="2021-12-16", version="1", precursor="true")
    root.append(l1_4)
    l1_4.text = sequence
    xmlstr = minidom.parseString(ET.tostring(root)).childNodes[0].toprettyxml(indent="    ")
    with open(file, 'a') as f:
        f.write(xmlstr)





with open(file, "w") as f:
    f.write('<?xml version="1.0" encoding="utf-8"?>\n')
    f.write('<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">\n')
for i in dict:
    accession = i.split("_")[1]
    fullname = i.split("_")[0]
    name = "".join([i[0] for i in fullname.split(" ")]).upper() + accession + "_SERMA"
    sequence = dict[i]
    GenerateXML(file, accession, name, fullname, sequence)
with open(file, "a") as f:
    f.write('<copyright> Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms Distributed under the Creative Commons Attribution (CC BY 4.0) License </copyright>\n')
    f.write('</uniprot>')

def BioUniprotIO_iterparser(PathUniprotXML: str, NS:str = "{http://uniprot.org/uniprot}"):
    """
    This function utilized Elementree iterparse to process Uniprot XML file in a way of each entry at a time
    Parameters
    ----------
    PathUniprotXML : str, directory of the XML file to be parsed
        DESCRIPTION. The default is PathUniprotXML.
    NS : str, fixed header for uniprot
        DESCRIPTION. The default is "{http://uniprot.org/uniprot}".
    Yields
    ------
    TYPE
        DESCRIPTION.
    """
    import xml.etree.ElementTree as etree
    from Bio.SeqIO.UniprotIO import Parser
    for event, elem in etree.iterparse(PathUniprotXML, events=("start", "end")):
        if event == "end" and elem.tag == NS + "entry":
            yield Parser(elem, return_raw_comments=False).parse()
            elem.clear()

def Write_fasta_from_XML(PathUniprotXML:str, PathFASTA:str, title = 'tr'):
    """
    This function is used to extract aa sequence from xml file and write into a fasta file
    ----------
    PathUniprotXML : str, directory of the XML file to be parsed
    PathFASTA : str, directory of the FASTA file to store the sequences
    ------
    return nothing
    """
    with open(PathFASTA, 'a') as fasta:
        for record in BioUniprotIO_iterparser(PathUniprotXML):
            seq = record.seq
            id = record.id
            name = record.name
            header = f">{title}|{id}|{name}\n"
            fasta.write(header + str(seq) + '\n')

Write_fasta_from_XML(PathUniprotXML='ATCC13880_hypothetical_proteins.xml', PathFASTA='test.fasta')
