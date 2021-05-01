# PTM predictor

Database search only retrieves existing information as hits. A poorly annotated database renders search results with low quality.

The PTM predictor is an application used to predict post-translational modification (PTM) site as candidates for poorly annotated proteins of bacteria based on sequence homology. These newly generated PTM site candidates will be added to existing database, which is formated as xml, and used for top-down mass spectral search. The goal of this application is to increase the number of hits at proteoform level during database searching. 

## Dependencies

PTM predictor relies on biopython and xml.etree.ElementTree to parse BLAST result and Uniprot xml file. It also utilizes matplotlib-venn to generate Venn diagram for comparison. These dependencies can be installed via [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/windows.html). The included environment.yml file makes this straightforward:
``` bash
conda env create -f environment.yml
conda activate PTM_predictor
```

## Running the application

Once we activate the virtual environment we can run the application by typing following codes in command line console with correct arguments
```bash
python PTM_predictor.py --query_xml \type\your\query xml\location\here --query_species [species name (e.g. Ecoli_K12)]  --sbjct_xml \type\your\subject xml\location\here --sbjct_species [species name (e.g. Ecoli_B)] --sl [length of short sequence]
```
**Note**: The elapsed time of the whole analysis depends on the size of the proteome. For the proteome that includes 5000 proteins, the elapsed time will be ~1 h


