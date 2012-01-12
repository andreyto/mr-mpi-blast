#!/bin/sh
/export/down/read-simul/mason/test/mason 454 -N 100 -o testSeq_454_400bp_mate_paired.fasta testSeq.fa -i -v -nm 400 -ne 0 -mp -rn 1
/export/down/read-simul/mason/test/mason  illumina -N 100 -o testSeq_illu_100bp_mate_paired.fasta testSeq.fa -i -n 100 -v -mp -rn 1

python /export/down/read-simul/make_query_for_mrblast.py testSeq_illu_100bp_mate_paired_1.fasta testSeq_illu_100bp_mate_paired_2.fasta testSeq_illu_100bp_mate_paired_mrblast.fasta
python /export/down/read-simul/make_query_for_mrblast.py testSeq_454_400bp_mate_paired_1.fasta testSeq_454_400bp_mate_paired_2.fasta testSeq_454_400bp_mate_paired_mrblast.fasta

sed 's/seq1/gi|9999|/g' testSeq_illu_100bp_mate_paired_mrblast.fasta > testSeq_illu_100bp_mate_paired_mrblast2.fasta
sed 's/seq2/gi|8888|/g' testSeq_illu_100bp_mate_paired_mrblast2.fasta > testSeq_illu_100bp_mate_paired_mrblast3.fasta
sed 's/seq1/gi|9999|/g' testSeq_454_400bp_mate_paired_mrblast.fasta > testSeq_454_400bp_mate_paired_mrblast2.fasta
sed 's/seq2/gi|8888|/g' testSeq_454_400bp_mate_paired_mrblast2.fasta > testSeq_454_400bp_mate_paired_mrblast3.fasta

python ../../../../src/tools/seqindexer/seqindexer.py testSeq_illu_100bp_mate_paired_mrblast3.fasta testSeq_illu_100bp_mate_paired_mrblast3.fasta.idx2
python ../../../../src/tools/seqindexer/seqindexer.py testSeq_454_400bp_mate_paired_mrblast3.fasta testSeq_454_400bp_mate_paired_mrblast3.fasta.idx2

#mpirun -np 4 ../../../../mrblast -db - -dbsize 12494903041 -evalue 1e-4 -num_threads 1 -window_size 0 -word_size 11 -searchsp 0 -num_descriptions 500 -num_alignments 10000   -penalty -5 -reward 4      -lcase_masking -dust yes -soft_masking true -export_search_strategy strategy_temp.txt  -max_target_seqs 2147483647
