"""
Created on Thur Apr 22 15:00:17 2021

@author: Taojunfeng Su
"""

#######################Functions#####################

import os

def get_PTMs_from_UniprotXML(PathUniprotAccessionNumber: str,
                             PathUniprotPTM: str,
                             PathUniprotXML: str,
                             Encoding="utf_8",
                             FeatureType=['modified residue', 'lipidation', 'glycosylation site'],
                             ):
    """
    Parameters
    ----------
    PathUniprotAccessionNumber : str, directory of a csv file containing Accession numbers and protein names
        DESCRIPTION. The default is PathUniprotAccessionNumber.
    PathUniprotPTM : str, directory of a csv file containing Accession numbers and PTMs info
        DESCRIPTION. The default is PathUniprotPTM.
    Encoding : TYPE, optional
        DESCRIPTION. The default is "utf_8".
    PathUniprotXML : str, directory of the XML file to be parsed
        DESCRIPTION. The default is PathUniprotXML.

    Returns
    -------
    None.

    """
    import codecs
    import csv
    with codecs.open(PathUniprotAccessionNumber, 'w', Encoding) as uniprot_accession_number, \
            codecs.open(PathUniprotPTM, "w", Encoding) as uniprot_PTM:
        AccessionNO_writer = csv.writer(uniprot_accession_number, quoting=csv.QUOTE_MINIMAL)
        PTM_writer = csv.writer(uniprot_PTM, quoting=csv.QUOTE_MINIMAL)

        #################write first line in each csv file################################
        AccessionNO_writer.writerow(['Accession_Number', 'Protein_Name'])
        PTM_writer.writerow(['Accession_Number', 'Feature_Key', 'PTM_description', 'Position', 'Length', 'Sequence'])

        #################process each entry (protein) in uniprot##########################
        for record in BioUniprotIO_iterparser(PathUniprotXML):
            AccessionNO_writer.writerow([record.id, record.name])
            for ft in record.features:
                if ft.type in FeatureType:
                    PTM_writer.writerow(
                        [record.id, ft.type, ft.qualifiers.get('description'),
                         list(ft.location)[0] + 1, len(record.seq), record.seq])


def BioUniprotIO_iterparser(PathUniprotXML: str, NS: str = "{http://uniprot.org/uniprot}"):
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


def NS_strip(elem_tag, NS: str = '{http://uniprot.org/uniprot}'):
    """
    This function is used to strip the uniprot print in a give tag
    :param tag:
    :param NS:
    :return:
    """
    return elem_tag.replace(NS, '')


def Extract_entry_in_xml_generatorStyle(PathUniprotXML: str, ProteinID: str):
    """
    This function is used to extract an entire entry for a given proteinID
    :param PathUniprotXML:
    :param ProteinID:
    :return: the whole entry that contains all information for a given protein. The returned data type is generator
    """
    import xml.etree.ElementTree as ET
    from Bio.SeqIO.UniprotIO import Parser
    for event, elem in ET.iterparse(PathUniprotXML, events=("start", "end")):
        if event == "end" and elem.tag == '{http://uniprot.org/uniprot}'+'entry' and elem.find('{http://uniprot.org/uniprot}'+'accession').text == ProteinID:
            yield Parser(elem, return_raw_comments=False).parse()
            elem.clear()


def Extract_entry_in_xml_returnStyle(PathUniprotXML: str, ProteinID: str):
    """
    This function is used to extract an entry for a given proteinID
    :param PathUniprotXML:
    :param ProteinID:
    :return: the whole entry that contains all information for a given protein. The returned data type is NOT generator
    """
    import xml.etree.ElementTree as ET
    for event, elem in ET.iterparse(PathUniprotXML, events=("start", "end")):
        elem.tag = NS_strip(elem.tag) # the uniprot address in each tag is removed
        if event == "end" and elem.tag == 'entry' and elem.find('accession').text == ProteinID:
            break
    return elem


def analyze_BLASTed_information_local_search(PathBLASTintermediate: str,
                                             Directory_intermediate_tobe_stored: str,
                                             E_value_threshold: float):
    """
    This function is used to analyze blasted information (xml format) generated from local search

    Parameters
    ----------
    E_value_threshold : TYPE, optional
        DESCRIPTION. The default is 0.001.

    Returns
    -------
    None.

    """
    from Bio.Blast import NCBIXML
    with open(PathBLASTintermediate) as results:
        result_handle = NCBIXML.parse(results)
        count = 1
        for blast_result in result_handle:
            # if there is no hit, we don't do anything and skip this iteration
            if len(blast_result.alignments) == 0:
                continue

            for alignment in blast_result.alignments:
                query_length = blast_result.query_length
                for hsp in alignment.hsps:
                    if hsp.expect < E_value_threshold:
                        Filename_blasted_information = "{protein_name}_{count}_blasted_result_info.txt".format(
                            protein_name=blast_result.query.split("|")[0], count=count)
                        PathBLASTinformation = os.path.join(Directory_intermediate_tobe_stored,
                                                            Filename_blasted_information)
                        count += 1
                        with open(PathBLASTinformation, 'w+') as blasted_fasta:
                            #######two conditions in PTM position calculation#####################
                            if int(blast_result.query.split("|")[2]) <= (query_length - 1) / 2:
                                Sbjct_PTM_position = int(
                                    blast_result.query.split("|")[2]) - hsp.query_start + hsp.sbjct_start
                            else:
                                Sbjct_PTM_position = int(
                                    (query_length - 1) / 2 - (hsp.query_start - 1) + hsp.sbjct_start)
                            Sbjct_PTM_aa = hsp.sbjct[Sbjct_PTM_position - hsp.sbjct_start]
                            blasted_fasta.write(
                                "*************************************Alignment**************************************" + '\n')
                            blasted_fasta.write("NCBI accession:" + alignment.accession + '\n')
                            blasted_fasta.write("Number of high scoring pair:" + str(len(alignment.hsps)) + '\n')
                            blasted_fasta.write("hit_id & definition:" + alignment.title + '\n')
                            blasted_fasta.write("length:" + str(alignment.length) + '\n')
                            blasted_fasta.write("e value:" + str(hsp.expect) + '\n')
                            blasted_fasta.write("score:" + str(
                                hsp.score) + '\n')  # hsp.bits for bit score;hsp.score for metric score
                            blasted_fasta.write("Sbjct_PTM_position:" + str(Sbjct_PTM_position) + '\n')
                            blasted_fasta.write("Sbjct_PTM_aa:" + Sbjct_PTM_aa + '\n')
                            blasted_fasta.write(
                                "Query" + "  " + "{0:<5}".format(str(hsp.query_start)) + "  " + hsp.query + "  " + str(
                                    hsp.query_start + hsp.align_length - 1) + '\n')
                            blasted_fasta.write(" " * 14 + hsp.match + '\n')
                            blasted_fasta.write(
                                "Sbjct" + "  " + "{0:<5}".format(str(hsp.sbjct_start)) + "  " + hsp.sbjct + "  " + str(
                                    hsp.sbjct_start + hsp.align_length - 1) + '\n')
                            blasted_fasta.write('\n')


def intermediate_result_summarize(PathBLASTintermediate: str,
                                  PathBLASTintermediateSummary: str,
                                  E_value_threshold: float):
    """
    The function is used to summarise information contained in BLAST intermediate xml file and stored in a csv file.
    The summary includes queried protein name, matched (subject) protein name, e-value, score, sbjct_PTM aa&position, Query_PTM aa&position
    :param PathBLASTintermediate: str
    :return: collected information
    """
    import csv
    from Bio.Blast import NCBIXML
    with open(PathBLASTintermediateSummary, 'w', newline="") as IntermediateStats, \
            open(PathBLASTintermediate) as results:
        Stats_writer = csv.writer(IntermediateStats, csv.QUOTE_MINIMAL)
        Stats_writer.writerow(
            ['Query_Accession_Number', 'Sbjct_Accession_Number', 'E_value', 'Match_Score', 'Query_sequence',
             'Sbjct_sequence', 'Query_PTM_aa', 'Query_PTM_position', 'Query_start', 'Query_end', 'Sbjct_PTM_aa',
             'Sbjct_PTM_position', 'Predicted_PTM', 'Query_Sbjct_aa_match'])
        result_handle = NCBIXML.parse(results)
        for blast_result in result_handle:
            # if there is no hit, we don't do anything and skip this iteration
            if len(blast_result.alignments) == 0:
                continue
            Query_Accession_Number = blast_result.query.split("|")[0]
            Predicted_PTM = blast_result.query.split("|")[1]

            for alignment in blast_result.alignments:
                query_length = blast_result.query_length
                for hsp in alignment.hsps:
                    if hsp.expect < E_value_threshold:
                        Sbjct_Accession_Number = alignment.title.split("|")[3]
                        E_value = hsp.expect
                        Match_Score = hsp.score
                        Query_sequence = hsp.query
                        Sbjct_sequence = hsp.sbjct
                        Query_start = hsp.query_start
                        Query_end = hsp.query_end
                        Query_PTM_position = blast_result.query.split("|")[2]
                        #######two conditions in PTM position calculation########################
                        if int(Query_PTM_position) <= (query_length - 1) / 2:
                            Query_PTM_aa = hsp.query[int(Query_PTM_position) - hsp.query_start]
                        else:
                            Query_PTM_aa = hsp.query[int((query_length - 1) / 2) - hsp.query_start + 1]

                        if int(blast_result.query.split("|")[2]) <= (query_length - 1) / 2:
                            Sbjct_PTM_position = int(
                                blast_result.query.split("|")[2]) - hsp.query_start + hsp.sbjct_start
                        else:
                            Sbjct_PTM_position = int((query_length - 1) / 2 - (hsp.query_start - 1) + hsp.sbjct_start)
                        Sbjct_PTM_aa = hsp.sbjct[Sbjct_PTM_position - hsp.sbjct_start]
                        #######evaluate if center aa in query sequence matches the one in subject sequence##########
                        if Query_PTM_aa == Sbjct_PTM_aa:
                            Query_Sbjct_aa_match = 'True'
                        else:
                            Query_Sbjct_aa_match = 'False'
                        Stats_writer.writerow(
                            [Query_Accession_Number, Sbjct_Accession_Number, E_value, Match_Score, Query_sequence,
                             Sbjct_sequence, Query_PTM_aa, Query_PTM_position, Query_start, Query_end, Sbjct_PTM_aa,
                             Sbjct_PTM_position, Predicted_PTM, Query_Sbjct_aa_match])


def local_BLAST(output_directory: str, Query_fasta: str, DB_tobe_searched_against: str,
                Query_species: str,
                Sbjct_species: str,
                Score_matrix_list=["PAM30"],
                E_value_threshold=0.0001):

    """
    This function is used to do local BLAST. The result will be an xml file containing all hits that passed the e value filter
    e.g. PAM30_Ecoli_MG1655_vs_Ecoli_B_BL21_DE3_intermediate.xml
    :param output_directory:
    :param Query_fasta:
    :param DB_tobe_searched_against:
    :param Query_species:
    :param Sbjct_species:
    :param Score_matrix_list:
    :param E_value_threshold:
    """
    import os
    from Bio.Blast.Applications import NcbiblastxCommandline
    for i in Score_matrix_list:
        Filename_blasted_intermediate = "{0}_{1}_vs_{2}_intermediate.xml".format(i, Query_species.replace(" ", "_"),
                                                                                 Sbjct_species.replace(" ", "_"))
        PathBLASTintermediate = os.path.join(output_directory, Filename_blasted_intermediate)
        blastp_cline = NcbiblastxCommandline(cmd='blastp', out=PathBLASTintermediate, matrix=i, outfmt=5,
                                             query=Query_fasta, db=DB_tobe_searched_against, evalue=E_value_threshold)
        blastp_cline()
    return Filename_blasted_intermediate, PathBLASTintermediate

def AccNO_PTM_position_extraction(Path_PTM_merge: str, Path_PTM_subject: str):
    """
    Description: This function is used to extract Accession Number and PTM position from both merged result and subject (strain that is used to search against).
    The result contains two lists of keys (tuple) that to be implemented in following dictionary.
    Format -- [(Sbjct_Accession_Number1, Sbjct_PTM_position1), (Sbjct_Accession_Number2, Sbjct_PTM_position2), .....].

    Warning: 1. The name of the column is expected to be fixed; 2. The Dataset_PTM_merge will be filtered based on column named Query_Sbjct_aa_match

    :param Dataset_PTM_merge:
    :param Dataset_PTM_subject:
    :return: one list of keys for each csv dataset
    """
    import pandas as pd
    Dataset_PTM_subject = pd.read_csv(Path_PTM_subject)
    Dataset_PTM_merge = pd.read_csv(Path_PTM_merge)

    # If the center aa in query is different from the one in sbjct, I will simply discard it
    Dataset_PTM_merge = Dataset_PTM_merge[Dataset_PTM_merge.Query_Sbjct_aa_match.eq(True)].reset_index(drop=True)

    merge_Acc_Position_key = []
    sbjct_Acc_Position_key = []
    for i in range(0, Dataset_PTM_merge.shape[0]):
        merge_Acc_Position_key.append(
            (Dataset_PTM_merge["Sbjct_Accession_Number"][i], Dataset_PTM_merge["Sbjct_PTM_position"][i]))
    for i in range(0, Dataset_PTM_subject.shape[0]):
        sbjct_Acc_Position_key.append((Dataset_PTM_subject["Accession_Number"][i], Dataset_PTM_subject["Position"][i]))
    return merge_Acc_Position_key, sbjct_Acc_Position_key


def Compare_list_of_Tuple(merge_Acc_Position_key: list, sbjct_Acc_Position_key: list):
    """
    Description: This function is used to compare two lists of tuple. Only the same tuple will be kept

    :param merge_Acc_Position_key:
    :param sbjct_Acc_Position_key:
    :return: a list of common Uniprot Acc NO. (Format -- [Acc1, Acc2, ....]) and their indexes in both dataset (Format -- (index in merge, index in sbjct))

    Note: the redundant elements were removed during the comparison
    """

    # The two lists may have redundant elements because same PTM site may have different modifications
    # So, the redundancy has to be removed before the comparison
    merge_Acc_Position_key = list(set(merge_Acc_Position_key))
    sbjct_Acc_Position_key = list(set(sbjct_Acc_Position_key))

    merge_sbjct_common = []
    merge_sbjct_pair_index = []
    for merge_i in range(0, len(merge_Acc_Position_key)):
        for sbjct_i, sbjct_v in enumerate(sbjct_Acc_Position_key):
            if sbjct_v == merge_Acc_Position_key[merge_i]:
                merge_sbjct_common.append(sbjct_v[0])  # the first element is Uniprot Acc NO.
                merge_sbjct_pair_index.append((merge_i, sbjct_i))

    return merge_sbjct_common, merge_sbjct_pair_index


def Make_two_circle_venn(NO_L: float, NO_R: float, NO_Inter: float, title: str, figsize=(8, 8), Char_L="Recorded", Char_R="Predicted", text_size=(22, 25, 30), gFile=True):
    from matplotlib_venn import venn2
    from matplotlib import pyplot as plt
    plt.figure(figsize=figsize)
    temp = venn2(subsets=(NO_L, NO_R, NO_Inter), set_labels=(Char_L, Char_R))
    for t in temp.set_labels:
        t.set_fontsize(text_size[0])
    for t in temp.subset_labels:
        t.set_fontsize(text_size[1])
    plt.title(title, fontsize=text_size[2])

    if gFile:
        plt.savefig(title+".png")
    else:
        pass


def Unique_element_in_tuple_list(primary_key_index: int, tuple_list: list):
    """
    Description: This function is used to grab unique elements from a tuple list based on one primary key
    :param primary_key_index:
    :param tuple_list:
    :return: a list with unique elements
    """
    try:
        type(tuple_list[0]) is tuple
    except:
        print("The list doesn't contain tuple")

    return list(set([x[primary_key_index] for x in tuple_list]))

def get_SeqInfo_from_UniprotXML(PathUniprotAccessionNumber: str,
                                ColumnOfAcc: int,
                                PathUniprotSeqInfo: str,
                                PathUniprotXML: str,
                                Encoding="utf_8"
                                ):
    """
    Parameters
    ----------
    PathUniprotAccessionNumber : str, directory of a csv file containing Accession numbers
        DESCRIPTION. The default is PathUniprotAccessionNumber.
    PathUniprotPTM : str, directory of a csv file containing Accession numbers and PTMs info
        DESCRIPTION. The default is PathUniprotPTM.
    Encoding : TYPE, optional
        DESCRIPTION. The default is "utf_8".
    PathUniprotXML : str, directory of the XML file to be parsed
        DESCRIPTION. The default is PathUniprotXML.

    Returns
    -------
    None.

    Notes :
     uniprot_accession_number file has a header

    """
    import codecs
    import csv
    with codecs.open(PathUniprotSeqInfo, "w", Encoding) as uniprot_SeqInfo, \
        codecs.open(PathUniprotAccessionNumber, "r", Encoding) as uniprot_accession_number:
        PathUniprotSeqInfo = csv.writer(uniprot_SeqInfo, quoting=csv.QUOTE_MINIMAL)
        PathUniprotAccessionNumber = csv.reader(uniprot_accession_number, quoting=csv.QUOTE_MINIMAL)
        #################write first line in each csv file################################
        PathUniprotSeqInfo.writerow(['Accession_Number', 'MW(Da)', 'Length', 'Sequence'])
        #################process each entry (protein) in uniprot##########################
        line_count = 0
        row_pre = '' # store previous Acc, if the next Acc is identical to previous one then just copy the previous info
        for row in PathUniprotAccessionNumber:
            if line_count == 0:  #skip the header
                line_count += 1
            else:
                if row_pre == row[ColumnOfAcc-1]:
                    PathUniprotSeqInfo.writerow(
                        [record.id, record.annotations.get('sequence_mass'), record.annotations.get('sequence_length'),
                         record.seq])
                else:
                    try:
                        for record in BioUniprotIO_iterparser(PathUniprotXML):
                            if row[ColumnOfAcc-1] == record.id:  #python starts with 0
                                PathUniprotSeqInfo.writerow([record.id, record.annotations.get('sequence_mass'), record.annotations.get('sequence_length'), record.seq])
                                row_pre = row[ColumnOfAcc - 1]
                                break
                    except:
                        print('Cannot extract info from XML')
def time_elapsed(func):
    import time
    def wrapper():
        print('Function is started...')
        before = time.time()
        func()
        print('Function is ended...' + '\n' + "Function took:", time.time()-before, "seconds")
    return wrapper
