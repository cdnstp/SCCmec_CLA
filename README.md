# Framework for Staphylococcal Cassette Chromosome *mec* Classification 

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

## Built With

* [prokka](add link) - Used to annotation
* [BLAST](add link) - Used to identified core elements
* [mash](add link) - Used to generate similarity networks

## Contributing

## Versioning

## Author

* **Felipe Sepúlveda** - *Initial work* - [fesepc](https://github.com/fesepc)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* **Ulab** [ugaldelab](https://github.com/ugaldelab) 

