#!/usr/bin/perl -w
## ==============================================================================
## Automated annotation tools
##
## program written by
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
##                                      http://weizhong-lab.ucsd.edu
## ==============================================================================

my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   $script_dir = "./" unless ($script_dir);

#### given scaffolds
#### ORFs are predicted by orf_finder from scaffolds
#### blastx query scaffolds against viral protein DB
#### this script select the ORFs confirmed by blastx

use Getopt::Std;
getopts("i:k:a:o:e:d:s:t:c:l:",\%opts);
die usage() unless ($opts{i} and $opts{d} and $opts{o});

my $ORF_file = $opts{i}; #### ORF by orf_finder
my $blx_file = $opts{d}; #### blastx alignment file in m8 format
my $output   = $opts{o}; #### output ORF
my $overlap_cutoff = $opts{c};
   $overlap_cutoff = 0.5 unless ($overlap_cutoff);
my $overlap_cutoff_l = $opts{l};
   $overlap_cutoff_l = 90 unless ($overlap_cutoff_l);

my ($i, $j, $k, $ll, $cmd);
my $output_looks_like = <<EOD; 

#query                                  subject          %       alnln   mis     gap     q_b     q_e     s_b     s_e     expect  bits
#0                                      1                2       3       4       5       6       7       8       9       10      11
sample|scaffold|164    NP_040516.2      99.16   1195    10      0       27104   30688   4       1198    0.0      2266   DNA polymerase [Human mastadenovirus C]
sample|scaffold|164    AP_000512.1      100.00  964     0       0       17038   14147   1       964     0.0      1936   hexon [Human adenovirus 1]
sample|scaffold|164    AP_000214.1      98.39   807     13      0       11781   9361    1       807     0.0      1380   100K [Human adenovirus 5]
sample|scaffold|164    NP_040518.2      97.98   645     12      1       25302   27233   4       648     0.0      1129   terminal protein precursor pTP [Human mastadenovirus C]
sample|scaffold|164    AP_000527.1      100.00  582     0       0       4798    3053    1       582     0.0      1084   fiber [Human adenovirus 1]
EOD

my %scaffold_blx_hits = ();
open(TMP, $blx_file) || die "can not open $blx_file";
while($ll = <TMP>){
  chop($ll);
  next if ($ll =~ /^#/);
  my @lls = split(/\s+/, $ll);

  my $sid   = $lls[0];
  my $seq_b = $lls[6];
  my $seq_e = $lls[7];
  my $frame_b = $seq_b; #### start of translation
  my $strand = "";

  if ($seq_b < $seq_e) {
    $strand = "+";
  }
  else {
    $strand = "-";
  }
  if (not defined($scaffold_blx_hits{$sid})) {
    $scaffold_blx_hits{$sid} = [];
  }
  ($seq_b, $seq_e) = sort {$a<=>$b} ($seq_b, $seq_e); ## sort in case of translated blast
  push(@{$scaffold_blx_hits{$sid}}, [$seq_b, $seq_e, $strand, $frame_b]);
}
close(TMP);


my $flag;
open(TMP, $ORF_file) || die "can not open $ORF_file";
open(OUT, "> $output") || die "can not write to $output";
while($ll = <TMP>){
  if ($ll =~ /^>/) {
#>ZYMOHIGH-VIRUS-MP-Rep_1|scaffold|164.19 /source=ZYMOHIGH-VIRUS-MP-Rep_1|scaffold|164 /start=7357 /end=7470 /frame=1 /length=38
    my ($ORF_id, $des) = split(/\s+/, substr($ll,1), 2);
    $ORF_id =~ s/\.(\d+)$/_$1/;
    my $sid = $ORF_id; $sid =~ s/_\d+$//;

    $flag = 0;
    my ($bb, $ee, $ff, $ss, $fb, $s1);
    if ($des =~ /start=(\d+)/) { $bb = $1;}
    if ($des =~ /end=(\d+)/)   { $ee = $1;}
    if ($des =~ /frame=(\S+)/) { $ff = $1; $ss = ($ff > 0) ? "+":"-"; }
    ($bb, $ee) = sort {$a<=>$b} ($bb, $ee);
    $fb = ($ff > 0) ? $bb : $ee;
    $s1 = ($ff > 0) ? 1   : -1;
    if (defined($scaffold_blx_hits{$sid})) {
      my @hits = @{$scaffold_blx_hits{$sid}};
      foreach $i (@hits) {
         my($seq_b, $seq_e, $strand, $frame_b) = @{$i};
         next unless ($ss eq $strand);
         next unless ( ($fb % 3) == ($frame_b % 3) );
         next unless ( overlap1($bb, $ee, $seq_b, $seq_e) > $overlap_cutoff * ($ee - $bb + 1));
         next unless ( overlap1($bb, $ee, $seq_b, $seq_e) > $overlap_cutoff_l );
         print STDERR "$ORF_id|begin|$bb|end|$ee|strand|$s1|fb|$fb\tmatch\t$sid|begin|$seq_b|end|$seq_e|strand|$strand|fb|$frame_b\n"; 
         $flag = 1;
         last;
      }
    }
    $ll = ">$ORF_id # $bb # $ee # $s1 $des" if ($flag);
  }
  print OUT $ll if ($flag);
}
close(TMP);
close(OUT);

sub overlap1 {
  my ($b1, $e1, $b2, $e2) = @_;
  return 0 if ($e2 < $b1);
  return 0 if ($b2 > $e1);
  return ( ($e1<$e2)? $e1:$e2 )-( ($b1>$b2)? $b1:$b2);
}


sub usage {
<<EOD;

The ORFs predicted by orf_finder are not all real genes,
This script confirms the real ORFs if the ORFs have good hit
to known protein database by blast.

Given a genomic DNA sequences, before running this script:
1) run orf_finder on the DNA to get raw ORFs
  orf_finder -i input_gDNA -o ORF-raw.faa -l 30 -L 30 -b 1 -e 1

2) run blastx on the DNA, against a known protein reference DB,
  blastx -query input_gDNA -out blastx_output -db prot_ref_db \
  -evalue 1e-6 -num_threads 4 -outfmt 6 -seg yes -max_target_seqs 5000

Then run this script as 

$script_name -i ORF-raw.faa -o ORF.faa -d blastx_output

  options:
    -i input ORF fasta file
    -o output ORF fasta file for blastx confirmed ORFs,
    -d blastx alignment file
    -c overlap_cutoff fraction of predicted ORF, default 0.5, 
       if the position of predicted ORF overlap with position of blastx alignment at this cutoff
       and if strand and read frame matches, this predicted ORF is confirmed
    -l overlap_cutoff length of overlap,  default 90 bp, 30 aa, 
       if the position of predicted ORF overlap with position of blastx alignment at this cutoff
       and if strand and read frame matches, this predicted ORF is confirmed
EOD
}
