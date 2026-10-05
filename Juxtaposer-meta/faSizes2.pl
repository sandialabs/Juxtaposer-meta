#! /usr/bin/perl -w
use strict;

my $usage = "perl $0 <FASTA file>";
die "usage:$usage\n" unless (@ARGV >0);
my ($genomefile) = @ARGV;
my $zip= 0 ;
my ($entry, $dna);
my $exit_status = system("gzip -t $genomefile 2> /dev/null");
if ($exit_status == 0) {
    $zip = 1;}

if ($zip) {open(IN, sprintf("zcat %s |", $genomefile))} else {open IN, $genomefile} 

while (<IN>) {
 chomp;
 if ( /^>(.+)/ ) {
  if ($entry) { Load() }
  $entry = $1;
  $dna = '';
 } else {$dna .= $_}
}
if ($entry) { Load() }

close(IN);

sub Load {
 $dna =~ s/[^a-zA-Z]//g;
 print "$entry\t", length($dna), "\n";
}



