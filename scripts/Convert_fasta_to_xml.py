'''
This script is used to convert fasta file to xml file
'''
import xml.etree.ElementTree as ET
import hashlib
import proteomic_tools as PT
from xml.dom import minidom
import math
from Bio import SeqIO
file = "test.xml"
dummy_index = 0
dict = {}
fasta_sequences = SeqIO.parse(open('LCIgSeqDB_noDecoys_nr.fasta'),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    dummy_index += 1
    dummy_acc = 'O' + str(dummy_index).zfill(5)
    dict[" ".join(name.split("_")[1:]).split('|')[0] + '_' + dummy_acc] = sequence

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
    name = "".join([i[0] for i in fullname.split(" ")]).upper() + accession + "_Human"
    sequence = dict[i]
    GenerateXML(file, accession, name, fullname, sequence)
with open(file, "a") as f:
    f.write('<copyright> Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms Distributed under the Creative Commons Attribution (CC BY 4.0) License </copyright>\n')
    f.write('</uniprot>')