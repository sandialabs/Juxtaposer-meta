#! /usr/bin/perl
use strict; use warnings;
use File::Basename;
use Getopt::Long;

my $usage = "Usage: perl $0 -g refPath -r readfile [options]\n";
die "$usage" unless @ARGV > 0;
my $help = "$usage\n".
	"-g <file> (Multi-) FastA DNA sequence file: path/filename\n".
	"-r <file> Path to interleaved FastQ reads\n".
	"-a <file> Path to existing bowtie2 alignment file (SAM or BAM) (default: None, bowtie is run freshly)\n".
	"-o <folder> Output directory (default: Current directory)\n".
	"-c <int> CPUs for programs with multithreading capabilities (>0) (default: 1)\n".
	"Additional options: -version, -usage, -help";

my $starttime= `date +%s.%N`;
chomp $starttime;
my $scriptdir = dirname(__FILE__);
my ($refPath, $readfile, $refpfx, $bowtie, $outDir, $cpu);
my ($palflag, $stdCommand) = 1, 0;
my $options_okay = GetOptions (
        'g=s' => \$refPath,
        'r=s' => \$readfile,
        'a=s' => \$bowtie,
	'o=s' => \$outDir,
	'c=i' => \$cpu,
        'usage' => sub { $stdCommand = $usage },
        'help' => sub { $stdCommand = $help } );

die "\n$stdCommand\n\n" if ($stdCommand);
die "\nUnknown arguments @ARGV\n\n$help\n\n" if (@ARGV);
die ("\n$help\n\n") unless $options_okay;
die "$usage" unless ( $refPath and $readfile);
die "Could not open file $refPath" unless (-e $refPath);
#die "Could not open file $readfile" unless (-e $readfile);
my @readarr=split(',', $readfile);
for (@readarr) {
	die "$_ does not exist\n" unless (-e $_)}

if ($outDir) { 
    unless ( -e $outDir ) { 
        unless ( system('mkdir -p', $outDir) == 0 ) {
            die "Could not make directory $outDir\n" 
        }
    }
    $outDir = `realpath $outDir`;
    chomp $outDir}
else {$outDir = `realpath .`;
    chomp $outDir}

die "$outDir is not a directory\n" unless ( -d $outDir );

if ($refPath =~ m{.*/([^/.]+)\..*$}) {
    $refpfx = "$outDir/$1";}
elsif ($refPath =~ m{([^/.]+)\..*$} ) {
    $refpfx = "$outDir/$1";}

$cpu = 1 unless $cpu;
die "cpu must be an integer\n" if $cpu =~ /[^\d]/;
my $posargpar = timepassed($starttime);

my $argtime=sprintf("%.3f",$posargpar);
print "TIMESTAMP: Arguments parsed in $argtime seconds\n";

my $curtime= `date +%s.%N`;
my $readlen = 0; # Maximum read length
if ( -e "${outDir}/read.lens") {
	open(IN,"${outDir}/read.lens");
	my $rl = <IN>;
	chomp $rl;
	my @line = split ",", $rl;
	if (scalar @line == 2 and -e $line[0] and `realpath $readarr[0]` eq `realpath $line[0]`) { 	
		$readlen = $line[1]}
	}
unless ( $readlen ) {
	my $zip =  0;
	my $exit_status = `gzip -t $readarr[0] 2>&1`;
	unless ($exit_status) {
		$zip = 1}
	# FIND MAXIMAL READLENGTH
	if ($zip) {open(IN, sprintf("zcat %s |", $readarr[0]))} else {open IN, $readarr[0]} 
	for (1..1000) {<IN>; my $seq = <IN>; last unless $seq; chomp $seq; my $len = length $seq; $readlen = $len if $len > $readlen; <IN>; <IN>}
	close IN;
	open(my $saverl,'>',"${outDir}/read.lens");
	print $saverl "$readarr[0],$readlen";
	close $saverl}
my $readsprep=timepassed($curtime);
print "TIMESTAMP: Readsize determined in $readsprep seconds\n";


# PROCESS REFERENCE GENOME
$curtime= `date +%s.%N`;
unless ( -e "${refpfx}.lens" 
	and -e "${refpfx}.mobile.bed" 
	and ( -e "{$bowtie}" or -e "${refpfx}.1.bt2") 
	and -e "${refpfx}.nin") {
		my $test_genome = `gzip -t $refPath 2>&1`; 
		unless ($test_genome) { # if any unfinished genome required tasks, check if compressed
			unless ( -e "${refpfx}.fa") { #only decompress file if it doesn't exist already
				warn "Genome is compressed, writing decompressed genome to ${refpfx}.fa\n"; 
				my $decompressed = system "zcat $refPath > ${refpfx}.fa";
				die "Could not decompress $refPath\n" unless ( $decompressed == 0)}
			warn "Using decompressed genome ${refpfx}.fa instead of $refPath\n";
			$refPath = "${refpfx}.fa"} }
my $genomeprep=timepassed($curtime);
print "TIMESTAMP: Genome decompression if compressed in $genomeprep seconds\n";

$curtime= `date +%s.%N`;
doCall("perl $scriptdir/faSizes.pl ${refPath} > ${refpfx}.lens") unless -f "${refpfx}.lens";
my $genomesize=timepassed($curtime);
print "TIMESTAMP: Genome sizes determined in $genomesize seconds\n";

$curtime= `date +%s.%N`;
my $mobile_cpu = $cpu -1;  # Hmmsearch with cpu==0 is safest
doCall("perl $scriptdir/find_mobile/mobilome_hunter.pl --fna ${refPath} --cpu $mobile_cpu") unless -f "$refpfx.mobile.bed";
my $mobhunt=timepassed($curtime);
print "TIMESTAMP: Mobilome annotations completed in $mobhunt seconds\n";

$curtime= `date +%s.%N`;
#SKIP IF USER SPECIFIES ALIGNMENT FILE
unless ($bowtie) {
	unless (-f "${refpfx}.1.bt2") {
		 doCall("bowtie2-build ${refPath} ${refpfx}")
		}
	doCall("bowtie2 --min-score L,-16.6,-0.073 -p $cpu --un ${outDir}/candidates.fq -x ${refpfx} -U ${readfile}* > /dev/null 2> ${outDir}/bowtie.results") unless -f "${outDir}/candidates.fq";
	doCall("perl $scriptdir/fq2fa.pl ${outDir}/candidates.fq ${outDir}/candidates.fa") unless -f "${outDir}/candidates.fa";
	unless (-f "${outDir}/ucandidates.results") {
		doCall("python $scriptdir/uniquereads.py ${outDir}/candidates.fa");
		system "mv ${outDir}/ucandidates.fa ${outDir}/candidates.fa"
		}
	}
#IF USER SPECIFIES ALIGNMENT FILE, PARSE OUTPUT FOR CANDIDATES ALIGNMENTS
else { 
	unless (-f "${outDir}/candidates.fa") {
		doCall("python $scriptdir/getcands.py ${bowtie} ${readfile} ${outDir} $readlen")
		}
	}
my $phase1=timepassed($curtime);
print "TIMESTAMP: Candidates reads collected from alignments in $phase1 seconds\n";
$curtime= `date +%s.%N`;
#PROCESS READS


doCall("makeblastdb -in ${refpfx}.fa -out ${refpfx} -dbtype nucl") unless -f "${refpfx}.nin";
doCall("blastn -num_threads $cpu -penalty -5 -perc_identity 95 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore btop' -query ${outDir}/candidates.fa -db ${refpfx} > ${outDir}/candidates.blast") unless -f "${outDir}/candidates.blast";
my $phase2=timepassed($curtime);
print "TIMESTAMP: BlastN run of candidate reads against genome completed in $phase2 seconds\n";
$curtime= `date +%s.%N`;
doCall("python $scriptdir/hitpairs.py ${outDir} ${refpfx} $readlen $palflag > ${outDir}/counts.txt") unless -f "${outDir}/juxtas.pkl";
doCall("python $scriptdir/network.py ${outDir}/juxtas.pkl $readlen");
my $hitp=timepassed($curtime);
print "TIMESTAMP: Final analysis completed in $hitp seconds\n";
my $totalruntime=timepassed($starttime);
print "TIMESTAMP: Total runtime was $totalruntime seconds\n";
warn "//FINISHED\n";
sub timepassed { 
 my $starttime = $_[0];
 $curtime = `date +%s.%N`;
 chomp $curtime;
 my $past = $curtime -$starttime;
 return $past;}

sub doCall {
 my $call = $_[0];
 warn "$call\n";  # Log the call
 my $code = system $call;
 die "Call failed, exiting: $call\nerror code $code\n" if $code;
}

