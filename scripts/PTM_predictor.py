"""
Created on Thurs Apr 22 18:53:50 2021

@author: Taojunfeng Su
"""
#!/usr/bin/python

import os
import sys
import shutil
import argparse
import pandas as pd
import xml.etree.ElementTree as ET
import PTM_predictor_module as PTMp
###########################define variables##################################################


@PTMp.time_elapsed
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--query_xml', type=str, default=r'~\..\uniprot-proteome_UP000000625.xml', help='Where is your query xml file?')
    parser.add_argument('--query_species', type=str, default='Ecoli_K12', help="What's the name of your query species (e.g., Ecoli_K12)")
    parser.add_argument('--sbjct_xml', type=str, default=r'~\..\uniprot-proteome_UP000002032.xml', help='Where is your subject xml file? The one you want to search against')
    parser.add_argument('--sbjct_species', type=str, default='Ecoli_B', help="What's the name of your sbjct species (e.g., Ecoli_B)")
    parser.add_argument('--sl', type=int, default=21, help='The length of each short sequence that tobe searched against your subject xml')
    args = parser.parse_args()
    #Current_path = os.path.dirname(args.query_xml)
    Current_path = os.getcwd()
    #e.g. Ecoli_MG1655
    Query_species = args.query_species
    #e.g. Ecoli_B_BL21_DE3
    Sbjct_species = DB_tobe_searched_against = args.sbjct_species

    #name of that uniprot xml that to be used as query---xml(query)
    #Filename_xml_query = os.path.basename(args.query_xml)  e.g. "proteome_UP000000625_Ecoli_k12_MG1655.xml"
    PathUniprotXML_query = args.query_xml
    #Filename_xml_sbjct = os.path.basename(args.sbjct_xml)   e.g. "proteome_UP000002032_Ecoli_B_BL21_DE3.xml"
    PathUniprotXML_sbjct = args.sbjct_xml

    search_length = args.sl
    ###########################define file name and directory##################################################

    #Accession number of each protein in xml(query)---Acc_csv(query)
    Filename_Accession_Number_csv_query = "Accession_Number_{}.csv".format(Query_species)
    PathUniprotAccessionNumber_query = os.path.join(Current_path, Filename_Accession_Number_csv_query)
    #PTM information of xml(query)---PTM_csv(query)
    Filename_PTM_info_csv_query = "PTM_information_{}.csv".format(Query_species)
    PathUniprotPTM_csv_query = os.path.join(Current_path, Filename_PTM_info_csv_query)
    #Updated xml(sbjct) with predicted PTM sites
    Filename_xml_sbjct_modified = 'uniprot-proteome_{}_PredicedPTM_Added.xml'.format(Sbjct_species)
    PathUniprotXML_xml_sbjct_modified = os.path.join(Current_path, Filename_xml_sbjct_modified)

    #Short sequences that generated based on PTM_csv(query)---short_sequence_fasta(query)
    Filename_short_sequence_fasta = "short_sequence_tobe_searched_{}.fasta".format(Query_species)
    PathShortFasta = os.path.join(Current_path, Filename_short_sequence_fasta)

    # The directory that results will be stored in. A folder will be made immediately
    Directory_intermediate_tobe_stored = os.path.join(Current_path, "BLAST_Intermediate_file\{}_intermediate".format(
        DB_tobe_searched_against))
    os.mkdir(os.path.join(Current_path, "BLAST_Intermediate_file"))
    os.mkdir(Directory_intermediate_tobe_stored)

    Filename_PTM_info_csv_sbjct = "PTM_information_{}.csv".format(Sbjct_species)
    PathUniprotPTM_csv_sbjct = os.path.join(Current_path, Filename_PTM_info_csv_sbjct)
    Filename_PTM_info_csv_predict = "PAM30_{0}_vs_{1}_intermediate_summary.csv".format(Query_species, Sbjct_species)
    PathUniprotPTM_csv_predict = os.path.join(Current_path, Filename_PTM_info_csv_predict)


    FeatureType = ['modified residue', 'lipidation', 'glycosylation site']

    Filename_Accession_Number_csv_sbjct = "Accession_Number_{}.csv".format(Sbjct_species)
    PathUniprotAccessionNumber_sbjct = os.path.join(Current_path, Filename_Accession_Number_csv_sbjct)

    ###########################run the main function##############################

    PTMp.get_PTMs_from_UniprotXML(PathUniprotAccessionNumber=PathUniprotAccessionNumber_query,
                             PathUniprotPTM=PathUniprotPTM_csv_query,
                             PathUniprotXML=PathUniprotXML_query)

    PTMp.get_PTMs_from_UniprotXML(PathUniprotAccessionNumber=PathUniprotAccessionNumber_sbjct,
                             PathUniprotPTM=PathUniprotPTM_csv_sbjct,
                             PathUniprotXML=PathUniprotXML_sbjct)

    with open(PathShortFasta, "w") as short_sequence_fasta:
        PTM_csv = pd.read_csv(PathUniprotPTM_csv_query)
        if search_length%2 == 0:
            print('search_length has to be an odd number')
            sys.exit()
        for i in range(0, PTM_csv.shape[0]):
            position = PTM_csv['Position'][i]
            length = PTM_csv['Length'][i]
            sequence = PTM_csv['Sequence'][i]
            description = PTM_csv['PTM_description'][i]

            ##determine how to cut the full length sequence
            # If the PTM position is smaller than the half of the search_length (aka close to N terminus),
            # the query sequence will be from the 1st aa to the 21st aa
            # e.g. |M|ESESDFAED....
            if position <= (search_length-1)/2:
                aa_end = search_length
                aa_start = 1
            # If the PTM position is close to C terminus,
            # the query sequence will end at the last aa and start from the aa that is 19 aa away from the end aa (if search_length = 21)
            # e.g. .....ARJF...|P|FRA
            elif length - position <= (search_length-1)/2:
                aa_end = length
                aa_start = aa_end - search_length + 1
            # the exception is that PTM position is not close to neither end
            else:
                aa_end = position + (search_length-1)/2
                aa_start = position - (search_length-1)/2
            short_seq_tobe_searched = sequence[int((aa_start-1)):int(aa_end)]
            accession = PTM_csv['Accession_Number'][i]
            #the header of the fasta(query) is formated as
            # >accession|PTM information|PTM position
            short_sequence_fasta.write(">" + accession + "|" + description + "|" + str(position) + "\n" + short_seq_tobe_searched + "\n")


    Filename_blasted_intermediate, PathBLASTintermediate = PTMp.local_BLAST(output_directory=Directory_intermediate_tobe_stored,
                                                                       Query_fasta=PathShortFasta, DB_tobe_searched_against=DB_tobe_searched_against,
                                                                       Query_species=Query_species, Sbjct_species=Sbjct_species,
                                                                       Score_matrix_list=["PAM30"], E_value_threshold=0.0001)


    # The summary of the output of BLAST---BLAST_intermediate_csv(all_sbjct)
    # The example of the header and first row are showing below
    # Query_Accession_Number	Sbjct_Accession_Number	E_value	Match_Score	Query_sequence	     Sbjct_sequence	   Query_PTM_aa	Query_PTM_position	Query_start	Query_end	Sbjct_PTM_aa	Sbjct_PTM_position	Predicted_PTM	        Query_Sbjct_aa_match
    #  P00509	                   A0A140ND68	        8.94E-17	154	  KELIVASSYSKNFGLYNERVG	KELIVASSYSKNFGLYNERVG	K	           246	             1	         21	         K	               246	         N6-(pyridoxal phosphate)lysine	    TRUE
    Filename_intermediate_summary = Filename_blasted_intermediate.replace('.xml', '_summary.csv')

    # concatenate directory with file name
    PathBLASTintermediateSummary = os.path.join(Directory_intermediate_tobe_stored, Filename_intermediate_summary)


    PTMp.analyze_BLASTed_information_local_search(PathBLASTintermediate=PathBLASTintermediate,
                                             Directory_intermediate_tobe_stored=Directory_intermediate_tobe_stored,
                                             E_value_threshold=0.0001)

    PTMp.intermediate_result_summarize(PathBLASTintermediate=PathBLASTintermediate,
                                  PathBLASTintermediateSummary=PathBLASTintermediateSummary,
                                  E_value_threshold=0.0001)


    #move the summarized file to current path for subsequent analysis
    shutil.copy(Directory_intermediate_tobe_stored+'\\' + Filename_intermediate_summary, Current_path+'\\'+Filename_intermediate_summary)

    try:
        #title of figures
        title_PTMsites = Sbjct_species + " PTM sites (PAM30)"
        title_Proteins = Sbjct_species + " Proteins with PTMs (PAM30)"

        # Making Venn diagram based on PTM site
        merge_Acc_Position_key, sbjct_Acc_Position_key = PTMp.AccNO_PTM_position_extraction(PathUniprotPTM_csv_predict, PathUniprotPTM_csv_sbjct)
        merge_sbjct_common, merge_sbjct_pair_index = PTMp.Compare_list_of_Tuple(merge_Acc_Position_key, sbjct_Acc_Position_key)
        PTMp.Make_two_circle_venn(NO_L=len(set(sbjct_Acc_Position_key))-len(merge_sbjct_pair_index), NO_R=len(set(merge_Acc_Position_key))-len(merge_sbjct_pair_index),
                             NO_Inter=len(merge_sbjct_pair_index), title=title_PTMsites, figsize=(13, 13))
    
        # Making Venn diagram based on protein
        unique_protein_names_sbjct = PTMp.Unique_element_in_tuple_list(0, sbjct_Acc_Position_key)  # some proteins may have multiple PTM sites. Here I only consider uniqueness at protein level
        unique_protein_names_merge = PTMp.Unique_element_in_tuple_list(0, merge_Acc_Position_key)
        unique_protein_names_inter = set(unique_protein_names_sbjct).intersection(unique_protein_names_merge)
    
        PTMp.Make_two_circle_venn(NO_L=len(unique_protein_names_sbjct)-len(unique_protein_names_inter), NO_R=len(unique_protein_names_merge)-len(unique_protein_names_inter),
                             NO_Inter=len(unique_protein_names_inter), title=title_Proteins, figsize=(13, 13))

    except:
        print("Seems like you don't have newly annotated PTMs")

    intermediate_csv_file = pd.read_csv(PathUniprotPTM_csv_predict, dtype={'Query_PTM_position': int, 'Query_start': int, 'Query_end': int, 'Sbjct_PTM_position': int})

    intermediate_csv_file_filtered = intermediate_csv_file.drop_duplicates(subset=['Sbjct_Accession_Number', 'Sbjct_PTM_position', 'Predicted_PTM'], keep='first').reset_index(drop=True)

    # simplify the information in the intermediate csv file by putting unique entries to a dictionary
    Dict_Acc_PTMpos_PTMdes = {} # The structure of this dictionary is {Acc1:{PTMpos1:[PTMdes1, PTMdes2], PTMpos2:[PTMdes2, PTMdes3]}, Acc2:{PTMpos3:[PTMdes1, PTMdes2]}}
    for i in range(0, intermediate_csv_file_filtered.shape[0]):
        Protein_Acc = intermediate_csv_file_filtered["Sbjct_Accession_Number"][i]
        PTM_position = int(intermediate_csv_file_filtered["Sbjct_PTM_position"][i])
        PTM_description = intermediate_csv_file_filtered["Predicted_PTM"][i]
        if Protein_Acc not in Dict_Acc_PTMpos_PTMdes:
            Dict_Acc_PTMpos_PTMdes[Protein_Acc] = {}
            Dict_Acc_PTMpos_PTMdes[Protein_Acc][PTM_position] = []
            Dict_Acc_PTMpos_PTMdes[Protein_Acc][PTM_position].append(PTM_description)
        elif Protein_Acc in Dict_Acc_PTMpos_PTMdes and PTM_position in Dict_Acc_PTMpos_PTMdes[Protein_Acc]:
            Dict_Acc_PTMpos_PTMdes[Protein_Acc][PTM_position].append(PTM_description)
        elif Protein_Acc in Dict_Acc_PTMpos_PTMdes and PTM_position not in Dict_Acc_PTMpos_PTMdes[Protein_Acc]:
            Dict_Acc_PTMpos_PTMdes[Protein_Acc][PTM_position] = []
            Dict_Acc_PTMpos_PTMdes[Protein_Acc][PTM_position].append(PTM_description)
        else:
            pass

    AccNO_csv_file = pd.read_csv(PathUniprotAccessionNumber_sbjct)
    with open(PathUniprotXML_xml_sbjct_modified, 'a') as f:
        # The is the main function that inserts new PTM sites into corresponding protein entry
        # logic:
        # start looping through each protein in subject species (e.g. Ecoli B BL21 DE3)
        # if the picked protein is not found in the Dict_Acc_PTMpos_PTMdes dictionary, then we just store that entry in the new xml file;
        # if the picked protein is found in the dictionary, we then test if the predicted PTM sites exist in the subject xml
        #  if they already exist, we will keep the entry as it is and store it in the newly written xml file;
        #  if any one of predicted PTM sites is new, we then create a new element for each new PTM sites, insert it/them into this entry, and store that modified entry in the newly written xml file.
        #  if there's the third possibility (I don't know yet), we will print the Acc of that protein for a manual check.
        
        #Add a header for ProsightPD database search
        f.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
        f.write('<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">\n')
  
        PathUniprotXML_sbjctAcceessionEntryDictionary = PTMp.Extract_access_entry_in_xml_dictionary(PathUniprotXML_sbjct)

        for i in range(0, AccNO_csv_file.shape[0]):
            Protein_Acc = AccNO_csv_file['Accession_Number'][i]
            entry = PathUniprotXML_sbjctAcceessionEntryDictionary[Protein_Acc]  # every protein will be extracted from the original xml no matter it has new PTM or not
            PTMpos = []
            if Protein_Acc in Dict_Acc_PTMpos_PTMdes.keys():
                for ft in entry.findall('feature'):
                    if ft.attrib['type'] in FeatureType:
                        PTMpos.append(int(ft.findall('location')[0].find('position').attrib['position']))

                # remove redundant PTM position for comparison
                PTMpos_original = set(PTMpos)
                PTMpos_predicted = set(Dict_Acc_PTMpos_PTMdes[Protein_Acc].keys())
                if not PTMpos or not PTMpos_predicted.issubset(PTMpos_original):
                    new_PTMpos_set = PTMpos_original.union(PTMpos_predicted) - PTMpos_original
                    for new_PTMpos in new_PTMpos_set:
                        for new_description in Dict_Acc_PTMpos_PTMdes[Protein_Acc][new_PTMpos]:
                            new_feature = ET.Element("feature", type='modified residue', description=new_description,
                                                     evidence='None', predicted="True")  # The type is fixed as modified residue, should be changed to broader range in the future
                            new_feature.text = '\n'
                            new_feature.tail = '\n'

                            new_sub1 = ET.SubElement(new_feature, "location")
                            new_sub1.text = '\n'
                            new_sub1.tail = '\n'

                            new_sub2 = ET.SubElement(new_sub1, "position", position=str(new_PTMpos))
                            new_sub2.text = '\n'
                            new_sub2.tail = '\n'

                            entry.insert(4, new_feature)

                    f.write(ET.tostring(entry, encoding="unicode"))
                elif PTMpos_predicted.issubset(PTMpos_original):
                    f.write(ET.tostring(entry, encoding="unicode"))  # No need to change anything, write as original
                else:
                    print(PTMpos_original)
                    print(PTMpos_predicted)
            else:
                f.write(ET.tostring(entry, encoding="unicode"))  # No need to change anything, write as original

        #Add several lines at the end for ProsightPD database search
        f.writelines('<copyright> Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms Distributed under the Creative Commons Attribution (CC BY 4.0) License </copyright>\n')
        f.writelines('</uniprot>')
if __name__ == '__main__':
    main()
