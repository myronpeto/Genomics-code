#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use IO::Zlib;
use Scalar::Util qw(looks_like_number);

###############################################################
#
#  Compares two sam files aligned from the same set of fastq
#  files and counts how many reads are aligned in different
#  locations. Used for comparing different versions of
#  alignment programs (bwa mem).
# 
#Author: Myron Peto 
my $version = "1.0.0";

#Last revision: February 19th - 2016
#
my $usage = "perl bam_compare.pl sam_file_1 sam_file_2 sam_outfile";
################################################################

sub usage(){
	print STDERR "$usage\n";
	exit;
}



# Three arguments passed to the script. The first two are the
# sam files to compare and the last is the sam file to output
# reads that disagree between the two alignments.
my $sam1 =  $ARGV[0];
my $sam2 =  $ARGV[1];
my $sam3 =  $ARGV[2];

# Temporary place holder for lines read from sam files
my $line1;

# Counts for matching/mismatching reads
my $count = 0;
my $tot = 0;
my $no_find = 0;
my $different = 0;


my $tmpid;

# Hash to hold the chromosome and position for each read id from
# the first sam file
my %readid = ();

my @parts1;
my @record1;


open (FILE1, "<$sam1") or die ("Can't upen the file $sam1");
print "The file opened is $ARGV[0]\n";

open (FILE2, "<$sam2") or die ("Can't upen the file $sam2");
print "The file opened is $ARGV[1]\n";

# The first do loop reads in reads from the first sam file
# and records their chromosome and start position in a hash.
# The key for the hash is the read id appended with a 1 or 2
# depending on whether it's read 1 or 2.

do {
   	$line1 = <FILE1>;
   	chomp($line1);
   	@parts1 = split /	/, $line1;
   	
   	# Test if it's read one
   	if ($parts1[1] & 64) {
   		$tmpid = $parts1[0] . "_1";
   		
	# Test if it's read two
   	} elsif ($parts1[1] & 128) {
   		$tmpid = $parts1[0] . "_2"
   	
   	# Otherwise something is wrong with the read	
   	} elsif ($line1) {
   		print "Invalid sam entry\n";
   	}
   	
   	# Record the information in a hash using the readid as the key
   	$readid{$tmpid} = "$parts1[3]	$parts1[1]";
   	$count++;
   	
   	# Print out progress every 20 million lines read
   	if ($count % 20000000 == 0) {
   		print $count, " reads read from file\n";
   	}
} until (!($line1));

# Open up for output the sam file to hold disagreeing reads.
open (OUT1, ">$sam3") or die("Can't open the file $sam3");
close (FILE1);

$count = 0;
$different = 0;
$tot = 0;
$no_find = 0;


do {
	$line1 = <FILE2>;
	if($line1) {
		chomp($line1);
		@parts1 = split /	/, $line1;
		
		# Check to see if it's read 1
		if ($parts1[1] & 64) {
	   		$tmpid = $parts1[0] . "_1";
	   		
	   	# Check to see if it's read 2
	   	} elsif ($parts1[1] & 128) {
	   		$tmpid = $parts1[0] . "_2";
	   		
	   	# Otherwise there's a problem with the read.
	   	} else {
	   		print "Invalid sam entry\n";
	   		print $parts1[1], "\n";
	   	}
	   	
	   	# Test to see if the read id was recorded from the first sam file
		if ($readid{$tmpid}) {
			@record1 = split /	/, $readid{$tmpid};
			
			# Test if both are read 1 and have the same start location.
			# If so, increment the count of matching reads.
			if ($record1[0] == $parts1[3] && ($record1[1] & 64) && ($parts1[1] & 64)) {
#				print OUT1 $line1, "\n";
				$count++;
				
			# Test if both are read 2 and have the same start location.
			# If so, increment the count of matching reads.
			} elsif ($record1[0] == $parts1[3] && ($record1[1] & 128) && ($parts1[1] & 128)) {
#				print OUT2 $line1, "\n";
				$count++;
				
			# Test if both are read 1 and have different start locations.
			# Print out the offending read and increment the count of non-matching reads.
			} elsif ($record1[0] != $parts1[3] && ($record1[1] & 64) && ($parts1[1] & 64)) {
				print OUT1 $line1, "\n";
				$different++;
				
			# Test if both are read 2 and have different start locations
			# Print out the offending read and increment the count of non-matching reads.
			} elsif ($record1[0] != $parts1[3] && ($record1[1] & 128) && ($parts1[1] & 128)) {
				print OUT1 $line1, "\n";
				$different++;
				
			# If none of the above finds a match, then the read must not have mapped in the
			# second alignment.
			} else {
				$no_find++;
			}
		}
		
		# Keep track of the total number of reads checked.
		$tot++;
		
		# Print out progress every 20 million lines read.
		if ($tot % 20000000 == 1) {
  			print $tot, " reads processed\t";
  			print "Didn't find: ", $no_find, "\t";
			print "Different location ", $different, "\t";
			print "Same location ", $count, "\n";
  		}
  	}
} until (!($line1));

# Print out totals for matching/mismatching reads.
print $tot, " reads processed\t";
print "Didn't find: ", $no_find, "\t";
print "Different location ", $different, "\t";
print "Same location ", $count, "\n";


close (OUT1);
close (FILE2);
    																				

