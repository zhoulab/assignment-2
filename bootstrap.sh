#!/bin/sh

module load python
virtualenv ve
. ve/bin/activate

pip install -r requirements.txt

wget http://geneontology.org/ontology/go-basic.obo -O ../data/
