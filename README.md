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

The basic logic of running the application is that **Run_PTM_predictor.bat** retrieves user defined parameters from **paramter.csv** file and passes them to the python scrpit to compile. This procedure can be achieved simply by runing the following code after we activate the virtual environment:
```bash
Run_PTM_predictor.bat
```



