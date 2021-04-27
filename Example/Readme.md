# Prepare inputs
## 1. Downloading Uniprot proteome XML file
In this example, we will use [_E. coli_ K12](https://www.uniprot.org/uniprot/?query=proteome:UP000000625) strain as a query and [_E. coli_ B](https://www.uniprot.org/uniprot/?query=proteome:UP000002032) strain as a subject in which predicted PTM sites will be added.

**Note**: to download xml file, click Download button, select format as XML, and hit Go button



## 2. Making a database for local BLAST search

BLAST+ is required for local BLAST search. Download [installer](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/) (.exe file) and install it in a proper directory. The example here uses version 2.8.1. Other versions may function as well but the BLAST result might be slightly different.
FASTA file of the subject species is also required for database making. Downloading the [FASTA](https://www.uniprot.org/uniprot/?query=proteome:UP000002032) file of _E. coli_ B strain and copy the file to the directory where you install BALST+

**Note**: to download fasta file, click Download button, select format as FASTA(canonical), and hit Go button

Open Command Prompt in the directory. Make sure you open it as administrator if it is installed in **C** drive. 
```bash
makeblastdb -in uniprot-proteome_UP000002032.fasta -dbtype prot -out E_coli_B_strain
```
Copy the three newly generated file to the folder where you want to store other inputs



## 3. Setting parameters for main function
The parameters will be stored in parameter.csv file, which is structured as following:
                  Current_path                 | Query_species | Sbjct_species |             Filename_xml_query             |             Filename_xml_sbjct           |   search_length
---------------------------------------------- | ------------- | ------------- | ------------------------------------------ | ---------------------------------------- | -------------
Directory where you want to store other inputs |    Ecoli_K12  |	   Ecoli_B   | uniprot-proteome_UP000000625_Ecoli_K12.xml | uniprot-proteome_UP000002032_Ecoli_B.xml | 21

**Note**: the order of the column in the csv file can't be changed, otherwise the script won't work
