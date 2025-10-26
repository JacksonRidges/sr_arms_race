#!/bin/zsh
#script that accepts beagle output file from ANGSD and runs PCANGSD
#-------------------------------------------------------
directory="/path/to/beagle"
beagle_name=""
output_name=""

cd ${directory}
pcangsd -b ${beagle_name} -o ${output_name} -t 8