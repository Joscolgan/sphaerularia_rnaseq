#!/bin/env bash
##################################################################################################
##
## Author: Joe Colgan                                            name:create_kallisto_reference.sh
##
## Date:
## 
##
## Purpose:
## This script takes a user-defined input FASTA sequence and outputs an indexed reference file
## for pseudoalignment with kallisto.  
## This script takes two arguments from the command line:
## 1) The name to provide to the index.
## 2) The path to input cDNA sequences provided in FASTA format.
##
##################################################################################################

## Take input arguments from the command line:
reference=$1
input=$2
## Check arguments exist.
## If no arguments provided, print usage to console.
if [[ $# -eq 0 ]] ; then
    echo 'No arguments provided'
    echo 'Usage: ./create_kallisto_reference.sh index_name path_to_input_sequences'
    exit 1
fi

## Generate kallisto-index:
./kallisto index -i "$reference" input/reference/"$input"
