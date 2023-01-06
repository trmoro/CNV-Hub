#!/bin/bash

echo "SARS Engine CNV-Hub BED Batch Entrypoint"
echo "$1 $2 $3 $4"
python3 -c'import cnv_app;cnv_app.process_bed("'"$1"'","'"$2"'","'"$3"'","'"$4"'")'
