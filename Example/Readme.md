# Prepare inputs
## 1. Downloading Uniprot proteome XML file
In this example, we will use [_E. coli_ K12](https://www.uniprot.org/uniprot/?query=proteome:UP000000625) strain as a query and [_E. coli_ B](https://www.uniprot.org/uniprot/?query=proteome:UP000002032) strain as a subject in which predicted PTM sites will be added.

**Note**: to download xml file, click Download button, select format as XML, and hit Go button



## 2. Making a database for local BLAST search

BLAST+ is required for local BLAST search. Download [installer](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/) (.exe file) and install it in a proper directory. The example here uses version 2.8.1. Other versions may function as well but the BLAST result might be slightly different.
FASTA file of the subject species is also required for database making. Downloading the [FASTA](https://www.uniprot.org/uniprot/?query=proteome:UP000002032) file of _E. coli_ B strain and copy the file to the directory where you install BALST+

**Note**: to download fasta file, click Download button, select format as FASTA(canonical), and hit Go button

Direct to the directory where you install BALST+, open Command Prompt, and use the following code to generate a blast database. The first parameter is the fasta file of subject strain. The second parameter is a user defined name of blast database. **Make sure you open it as administrator if it is installed in C drive. Otherwise there won't be any output** 
```bash
makeblastdb -in uniprot-proteome_UP000002032.fasta -dbtype prot -out Ecoli_B
```
Copy the three newly generated file to the folder where you want to store other inputs


## 3. Setting arguments
As we described in front page, we need to set five arguments before runnning the application. In this example we can set arguments as following.
--query_xml ~\uniprot-proteome_UP000000625.xml
--query_species Ecoli_K12
--sbjct_xml ~\uniprot-proteome_UP000002032.xml
--sbjct_species Eoli_B
--sl: 21

**Note**: 
1. The nomenclature of the query_species is not rigorious but sbjct_species must match the name of BLAST index files (phr, pin and psq) you created. You can put any name as long as it's informative
2. The search length must be an odd number. This requirment enables one amino acid locates at the center of the sequence
