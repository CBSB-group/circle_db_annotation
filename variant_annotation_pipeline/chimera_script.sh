#!/bin/bash
export PATH="/home/damiano/Chimera/bin:$PATH"


# Controlla se sono stati forniti abbastanza argomenti
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <pdb_file> <chain> <position> <distance> <output_file>"
    exit 1
fi

# Estrai gli argomenti
pdb_file=$1
chain=$2
position=$3
distance=$4
output_file=$5

# Esegui i comandi in Chimera
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

