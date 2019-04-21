


## Define reference file:
reference=$1
input=$2

## Generate kallisto-index:
./kallisto index -i "$reference" input/reference/"$input"
