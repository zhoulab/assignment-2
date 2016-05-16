#!/bin/sh

module load python

source ve/bin/activate
pip freeze | xargs pip uninstall -r requirements.txt
rm -r ve/
rm -r go-basic.obo