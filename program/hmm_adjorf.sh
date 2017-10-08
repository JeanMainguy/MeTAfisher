#!/bin/bash

# echo 'START HMM'
hmmdb=dependence/ALL_plus_MET_curatted.hmm
seqfile=data/temporary_adjOrf.faa
out_table=./data/table_hmm.txt
# echo HMMSEARCH
hmmsearch -E 0.5 -o tmp.txt  --domtblout ${out_table}  ${hmmdb} ${seqfile}
# echo 'END HMM'
