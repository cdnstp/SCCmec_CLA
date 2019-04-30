# Pipeline for Staphylococcal Cassette Chromosome *mec* Classification 

## Prerequisites

* python 2.7 (2.7.12)
* virtualenv

## Installing

### Create working dir

```
$ mkdir workingspace
$ cd workingspace
```

### Create a virtualenv & activate it
```
$ virtualenv sccmec_cla
$ source ./sccmec_cla/bin/activate

```

### Download and unzip 
```
$ unzip SCCmec_CLA-master.zip
$ cd SCCmec_CLA-master

```
### Setup PROKKA database
```
$ ./helpers/programs/prokka-1.12/bin/prokka --setupdb

```
### Install requirements
```
$ pip install -r requirements.txt
```

## Usage
### Put sequences (contigs files) inside INPUT folder and run

```
$ python sccmec_classification.py
```

## OUTPUT

* attL_<CONTIG_FILENAME>.fasta <- left end sccmec
* attR_<CONTIG_FILENAME>.fasta <- right end sccmec
* sccmec_<CONTIG_FILENAME>.fasta <- cassette fasta format
* sccmec_<CONTIG_FILENAME>_type.txt <- current format annotation
* core_elements_sccmec_<CONTIG_FILENAME>.txt <- detail current format annotation
* annotation_table_sccmec_<CONTIG_FILENAME>.txt <- sccmec annotation table 
* sccmec_<CONTIG_FILENAME>_neighbors_cassettes.txt <- close related cassettes
* fig.png <- network viz using matplotlib
* sccmec_<CONTIG_FILENAME>_cytoscape_network.sif <- Cytoscape format files
* sccmec_<CONTIG_FILENAME>_cytoscape_network.eda <- Cytoscape format files
* SCCmec_<CONTIG_FILENAME>.png <- Graphical Representation 



## Built With

* [prokka](add link) - Used to make annotation
* [BLAST](add link) - Used to identify core elements
* [mash](add link) - Used to generate similarity networks
* [DnaFeaturesViewer](add link) - Used create a graphical representation

## Contributing

## Versioning

## Author

* **Felipe Sepúlveda** - *Initial work* - [fesepc](https://github.com/fesepc)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* **Ulab** [ugaldelab](https://github.com/ugaldelab) 

