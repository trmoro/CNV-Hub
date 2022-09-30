#!/bin/bash

echo "SARS CNV Hub"
cd XCNV
sh Install.sh
cd ..
python3 xcnv_mod.py

if [ $1 = "test_xcnv" ]
then 
	echo "Test XCNV"
	python3 -c'import cnv_batch;cnv_batch.test_xcnv()'
fi

if [ $1 = "process" ]
then 
	echo "Process BED File"
	python3 -c'import cnv_batch;cnv_batch.process("'"$2"'")'
fi
