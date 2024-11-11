# Prepare inputs
## 1. Downloading Uniprot proteome XML file
In this example, we will use [_E. coli_ K12](https://www.uniprot.org/uniprot/?query=proteome:UP000000625) strain as a query and [_E. coli_ B](https://www.uniprot.org/uniprot/?query=proteome:UP000002032) strain as a subject in which predicted PTM sites will be added.

**Note**: to download xml file, click Download button, select format as XML, and hit Go button



## 2. Making a database for local BLAST search
There are two methods to create BLAST databases for the BLAST search. Both method require installed [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/) and a fasta file. The fasta file of the subject species can be downloaded from UniProtKB. In this example, download [FASTA](https://www.uniprot.org/uniprot/?query=proteome:UP000002032) file of _E. coli_ B strain.

**Note**: to download a fasta file, click Download button, select format as FASTA(canonical), and hit Go button.
1. Use Python script via Commond Promt
   Type the following codes in your terminal.
   ```bash
   python MakeBLASTdb.py --fasta \type\your\fasta file\location\here --result_name [species name (e.g. Ecoli_K12)]
   ```
   
2. Use BLAST+. 
   Copy the fasta file of _E. coli_ B strain to the directory where you install BALST+.
   Direct to the directory where you install BALST+, open Command Prompt, and use the following code to generate a blast database. The first parameter is the fasta file of subject strain. The second parameter is a user defined name of blast database. **Make sure     you open it as administrator if it is installed in C drive. Otherwise there won't be any output** 
   ```bash
   makeblastdb -in uniprot-proteome_UP000002032.fasta -dbtype prot -out Ecoli_B
   ```
Both methods will generate three files, .pin, .psq, and .phr. Copy these files to the directory storing all the codes.


## 3. Setting arguments
As we described in front page, we need to set five arguments before runnning the application. In this example we can set arguments as following.
--query_xml ~\uniprot-proteome_UP000000625.xml
--query_species Ecoli_K12
--sbjct_xml ~\uniprot-proteome_UP000002032.xml
--sbjct_species Eoli_B
--sl: 21

**Note**: 
1. The nomenclature of the query_species is not rigorious but sbjct_species must match the name of BLAST index files (phr, pin and psq) you created.
2. The search length must be an odd number. This requirment enables one amino acid locates at the center of the sequence
