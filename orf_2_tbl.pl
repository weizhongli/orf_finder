#!/usr/bin/perl

print "#ORF\tsource_DNA\tstart\tend\tstrand\tlength\n";
while($ll=<>){
  next unless ($ll =~ /^>/);
  $ll =~ s/\s+$//g;
  $orf = $ll; $orf =~ s/\s+.+$//;
  $orf = substr($orf,1);

  my $source = "-";
  my $start = "-";
  my $end = "-";
  my $strand = "=";
  my $length = "-";
  if ($ll =~ /\/source=(\S+)/) {$source=$1;}
  if ($ll =~ /\/start=(\S+)/ ) {$start =$1;}
  if ($ll =~ /\/end=(\S+)/   ) {$end   =$1;}
  if ($ll =~ /\/length=(\S+)/) {$length=$1;}
  if ($ll =~ /\/frame=(\S+) /) {$strand=$1;  $strand=($strand>0)?"+":"-";}

  print "$orf\t$source\t$start\t$end\t$strand\t$length\n";
}
