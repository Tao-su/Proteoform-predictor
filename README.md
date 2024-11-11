# Proteoform-predictor

Database search only retrieves existing information as hits. A poorly annotated database renders search results with low quality.

The Proteoform-predictor is an application used to predict post-translational modification (PTM) site as candidates for poorly annotated proteins of bacteria based on sequence homology. These newly generated PTM site candidates will be added to existing database, which is formated as xml, and used for top-down mass spectral search. The goal of this application is to increase the number of hits at proteoform level during database searching. 

## Downloading scripts

Before using, users can download a zip file containing all codes and files into their local PCs. Make sure the directory storing these files is the working directory. 

## Installing Dependencies
Besides Python and its modules, Proteoform-predictor requires [BLAST+] (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/) to work. The example here uses version 2.8.1. Other versions may function as well but the BLAST result might be slightly different. In order to use functions related to BLAST+, be sure to add ~BLAST+\ncbi-blast-2.8.1\bin into the environment variable. The exact steps of adding environment variables can vary a lttile in different versions of Windows system, but it can be found online. 
Proteoform-predictor relies on biopython and xml.etree.ElementTree to parse BLAST result and Uniprot xml file. It also utilizes matplotlib-venn to generate Venn diagram for comparison. These dependencies can be installed via [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/windows.html) using the following codes in Command Prompt. Before execute the codes, make sure the working directory contains the environment.yml file.
``` bash
conda env create -f environment.yml
conda activate PTM_predictor
```
To properly create an environment, make sure add the ~Anaconda/ and ~Anaconda/Library/bin into the environment variable. If users encounter issues during the activation of the environment, try the following code.
``` bash
conda init cmd.exe
```

## Running the application

Once we activate the virtual environment, we can run the application by typing following codes Command Prompt with correct arguments. The instruction of downloading required XML files and using BLAST+ can be found in the Example folder. Once the required files are prepared, users can use the following cmd line to execute the program. In brief, --query_xml and --query_species take the file location of the query xml and the name of the query sepcies, respectively. --sbjct_xml and --sbjct_species also require a file location and species name but they represent the subject database. Finally, the search length represents the length of each short sequence. The default length is 21.

```bash
python PTM_predictor.py --query_xml \type\your\query xml\location\here --query_species [species name (e.g. Ecoli_K12)]  --sbjct_xml \type\your\subject xml\location\here --sbjct_species [species name (e.g. Ecoli_B)] --sl [length of short sequence]
```
To check what arguments should be given to the program, we can call help function as follow

```bash
python PTM_predictor.py -h
```

**Note**: The elapsed time of the whole analysis depends on the size of the proteoform database. For the E. coli K12 proteome, the elapsed time will be ~1 h


