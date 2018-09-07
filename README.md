# orf_finder

This is not a gene prediction program. It just calls ORFs from

genomic sequences. The output ORFs may not be real genes. So it 

is intended to use together with other method, such as blast against

a reference database to validate the ORFs.


# To compile:

make clean

make

To run orf_finder:

./orf_finder -h for help

./orf_finder -l 30 -L 30 -t 11 -i test1.fna -o test1.faa

./orf_2_tbl.pl < test1.faa > test1.tab

