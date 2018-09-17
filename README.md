# orf_finder
This is not a gene prediction program. It just calls ORFs from
genomic sequences. The output ORFs may not be real genes. So it 
is intended to use together with other method, such as blast against
a reference database to validate the ORFs.


# To compile:
```
make clean
make
```

# To run orf_finder:
```
./orf_finder -h for help
./orf_finder -l 30 -L 30 -t 11 -i test1.fna -o test1.faa
./orf_2_tbl.pl < test1.faa > test1.tab
```

# Optional, further validation based on blast alignment to known genes
The ORFs predicted by orf_finder are not all real genes,
script select_ORF_by_blastx.pl confirms the real ORFs if the ORFs have good hit
to known protein database by blast.

Given a genomic DNA sequences, before running this script:
1. run orf_finder on the DNA to get raw ORFs
```
orf_finder -i input_gDNA -o ORF-raw.faa -l 30 -L 30 -b 1 -e 1
```
2. run blastx on the DNA, against a known protein reference DB,
```
blastx -query input_gDNA -out blastx_output -db prot_ref_db -evalue 1e-6 -num_threads 4 -outfmt 6 -seg yes -max_target_seqs 5000
```
3. then run this script as
```
select_ORF_by_blastx.pl -i ORF-raw.faa -o ORF.faa -d blastx_output
```
