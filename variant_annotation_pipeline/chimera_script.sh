#!/bin/bash

#check arguments

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <pdb_file> <chain> <position> <distance> <output_file>"
    exit 1
fi

# arguments extraction
pdb_file=$1
chain=$2
position=$3
distance=$4
output_file=$5

# chimera command

chimera --nogui << EOF
open $pdb_file
zone :$position.$chain $distance
select up
namesel surroundings
select :$position.$chain
namesel residue
select surroundings | residue
write selected 0 $output_file
stop
EOF

