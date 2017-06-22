# Framework for Staphylococcal Cassette Chromosome *mec* Classification 

## Prerequisites

python 2.7 (2.7.12)
virtualenv

## Installing

```
$ mkdir workingspace
$ cd workingspace
$ virtualenv sccmec_cla
$ source ./sccmec_cla/bin/activate
$ unzip SCCmec_CLA-master.zip
$ cd SCCmec_CLA-master
$ ./helpers/programs/prokka-1.12/bin/prokka --setupdb
$ pip install -r requirements.txt
```

## Usage

```
$ python sccmec_classification.py
```

## Built With

* [prokka](add link) - Used to annotation
* [BLAST](add link) - Used to find homologues
* [mash](add link) - Used to generate similarity networks

## Contributing

## Versioning

## Author

* **Felipe Sepúlveda** - *Initial work* - [fesepc](https://github.com/fesepc)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Ulab (add link)


