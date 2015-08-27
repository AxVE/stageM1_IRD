#!/usr/bin/perl
#Write by A. Verdier (axel.l.verdier@free.fr)
#
#v. 20150827

use strict;
use warnings;

sub help{
	print "This script is a pipeline to run multiple stats on alignment result of paired reads.\n";
	print "To use it : ./analyseAlignment.pl <file.bam> <file1.fq.gz> <file2.fq.gz>\n";
	print "Before that, set the correct values to the user constants\n";
	print "To display this help: ./analyseAlignment.pl -h\n";

	exit 0;
}


# User constants (must be edited for the project)
use constant {
	#Server informations
	LOGIN => "averdier",
	SERVER => "nas2",
	#File/folder
	FILE_REF => "/data/projects/evoreprice/Reference/Osativa_204_v7.0.fa",
	FILE_REF_INDEX => "/data/projects/evoreprice/Reference/Osativa_204_v7.0.fa.fai",
	FILE_ANNOT => "/data/projects/evoreprice/Reference/Osativa_204_v7.0.gene_exons.cuffcomp.rRNAremoved.gtf",
	#Script parameters
	RMV_ENV => 1, #Boolean if work files and dir must be remove (0 = false)
};

# Script constants
use constant {
	#Tools
	SAMTOOLS => "/usr/local/bin/samtools",
	#File/folder
	FILES => "Files/",
	DIR_OUTPUT => "out-analyseAlignment/",
	SCRIPT_GAPSTATS => "/data/projects/evoreprice/Scripts/AnalyseAlignment/gap_stats.pl",
	SCRIPT_COUNTMMBINNED => "/data/projects/evoreprice/Scripts/AnalyseAlignment/count_mm_binned.pl",
	SCRIPT_GETGAPMMPOS => "/data/projects/evoreprice/Scripts/AnalyseAlignment/get_gap_mm_pos.pl",
};

### MAIN ###

#prepare directories
print "Prepare work environment\n";
system("mkdir -p ".(DIR_OUTPUT)." ".(FILES));
system("mkdir ".(FILES)."Fastq/");
my $dir_gap_stats=(DIR_OUTPUT)."gapstats/";
system("mkdir $dir_gap_stats");

## parse and check parameters ##
print "Parse parameters\n";
if($ARGV[0] eq "-h"){help();} #Display help and stop

if(@ARGV != 3){
	print STDERR "Error: expected 3 parameter: <file.bam> <file1.fq.gz> <file2.fq.gz>\n";
	help(); #Display help and stop
}

(my $bamFile) = $ARGV[0]=~ m/^.*\/(.*)$/; #Get file name + ext

## get files into (FILES) folder##
print "Get files:\n";
#note : scp return 0 for OK
#copy reference and ref index
print("\t".(FILE_REF)."\n");
system("scp -q ".(LOGIN)."@".(SERVER).":".(FILE_REF)." ".(FILES)."ref.fa")==0  or die("Can't get ".(FILE_REF));
print("\t".(FILE_REF_INDEX)."\n");
system("scp -q ".(LOGIN)."@".(SERVER).":".(FILE_REF_INDEX)." ".(FILES)."ref.fa.fai")==0 or die("Can't get ".(FILE_REF_INDEX));
print("\t".(FILE_ANNOT)."\n");
system("scp -q ".(LOGIN)."@".(SERVER).":".(FILE_ANNOT)." ".(FILES)."annot.gtf")==0 or die("Can't get ".(FILE_ANNOT));
print("\t".(SCRIPT_GAPSTATS)."\n");
system("scp -q ".(LOGIN)."@".(SERVER).":".(SCRIPT_GAPSTATS)." ".(FILES)."gap_stats.pl")==0 or die("Can't get ".(SCRIPT_GAPSTATS));
print("\t".(SCRIPT_COUNTMMBINNED)."\n");
system("scp -q ".(LOGIN)."@".(SERVER).":".(SCRIPT_COUNTMMBINNED)." ".(FILES)."count_mm_binned.pl")==0 or die("Can't get ".(SCRIPT_COUNTMMBINNED));
print("\t".(SCRIPT_GETGAPMMPOS)."\n");
system("scp -q ".(LOGIN)."@".(SERVER).":".(SCRIPT_GETGAPMMPOS)." ".(FILES)."get_gap_mm_pos.pl")==0 or die("Can't get ".(SCRIPT_GETGAPMMPOS));
#copy bam file
print("\t$ARGV[0]\n");
system("scp -q ".(LOGIN)."@".(SERVER).":".$ARGV[0]." ".(FILES)."$bamFile")==0 or die("Can't get $ARGV[0]");
#copy fastq files
print("\t$ARGV[1]\n");
system("scp -q ".(LOGIN)."@".(SERVER).":".$ARGV[1]." ".(FILES)."Fastq/reads_1.fq.gz")==0 or die("Can't get $ARGV[1]");
print("\t$ARGV[2]\n");
system("scp -q ".(LOGIN)."@".(SERVER).":".$ARGV[2]." ".(FILES)."Fastq/reads_2.fq.gz")==0 or die("Can't get $ARGV[2]");

## generate md line in bam (neede by the analyse scripts) ##
my $bamName = $bamFile;
$bamName=~s/\.\w+$//; #Remove extension

print "Generate md line\n";
	#Store the generate bam in $bamMdFile as: bamFile=name.bam => bamMdFile=name_md.bam
my $bamMdFile = $bamFile;
$bamMdFile =~ s/^(.*)\.(.*)$/$1_md\.$2/;
	#cmd line: samtools calmd -b <in.bam> <ref.fa> [> out.bam]
my $cmd=(SAMTOOLS)." calmd -b ".(FILES).$bamFile." ".(FILES)."ref.fa > ".(FILES).$bamMdFile;
print "\t$cmd\n";
system($cmd)==0 or die ("Error $?");

## gap-stats ##
print "Script: gap_stats.pl\n";
	#Run script
$cmd=(SAMTOOLS)." view ".(FILES).$bamMdFile." | ".(FILES)."gap_stats.pl ".(FILES)."annot.gtf ".$bamName."_gapstats";
print "\t$cmd\n";
system($cmd)==0 or die ("Error $?");
	 #Move results files
print "\tMove results files\n";
system("mv os/ oc/ iac/ is/ $dir_gap_stats")==0 or die ("Failed to move results file of gap_stats.pl");

## count_mm_binned ##
print "Script: count_mm_binned.pl\n";
	#Run script
$cmd=(FILES)."count_mm_binned.pl ".(FILES)."Fastq/ ".(FILES).$bamMdFile." 33 > ".(DIR_OUTPUT).$bamName."_countmmbinned"; # 33 is the quality score based in sam/bam
print "\t$cmd\n";
system($cmd)==0 or die ("Error $?");

## get_gap_mm_pos ##
print "Script: get_gap_mm_pos.pl\n";
	#Run script
$cmd=(SAMTOOLS)." view ".(FILES).$bamMdFile." | ".(FILES)."get_gap_mm_pos.pl > ".(DIR_OUTPUT).$bamName."_getgapmmpos";
print "\t$cmd\n";
system($cmd)==0 or die ("Error $?");

## Remove work files and dir
if((RMV_ENV)){
	print "Remove files and directories\n";
	print "\tRemove Files\n";
	system("rm -r ".(FILES))==0 or die ("Problem during deleting of ".(FILES));
}
