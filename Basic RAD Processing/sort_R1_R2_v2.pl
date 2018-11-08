use strict;
use warnings;

# This script 
#	reads two fasta files and matches R1 and R2 reads for the same loci in order into two new files
#	R1 loci without R2 matches are reported


# run: perl sort_R1_R2.pl <R1_fasta_file> <R2_fasta_file>

###################################


my $R1File = "$ARGV[0]";
my $R2File = "$ARGV[1]";

#read in R1 file and put loci into an array
my @master;
open IN1, "<$R1File";
foreach (<IN1>) {
	chomp $_;
	if ($_) { #if the line is not blank
		if ($_ =~ m/^>/) { 
			my @parts = split / /, $_; #split off end of header for storage to master
			push @master, $parts[0]; 
		} 
	}
}
close IN1;


#read in R2
my %hash;
my $header;
open IN2, "<$R2File";
foreach (<IN2>) {
	chomp $_;
	if ($_) { #if the line is not blank
		if ($_ =~ m/^>/) { $header =  $_; }
		elsif ($_) { 
			my @parts = split / /, $header;
			$hash{$parts[0]} = "$parts[1]\n$_"; #collect header base, and then its end with sequence information
		} 
	}
}
close IN2;


#print R2
my %missing;
open OUT2, ">$R2File.sorted";
foreach my $entry (@master) { 
	if ($hash{$entry}) { print OUT2 "$entry $hash{$entry}\n"; }
	else { $missing{$entry}++; } 
}
close OUT2;


#read R1 again
%hash = ();
open IN1, "<$R1File";
foreach (<IN1>) {
	chomp $_;
	if ($_) { #if the line is not blank
		if ($_ =~ m/^>/) { $header = $_; }
		elsif ($_) { 
			my @parts = split / /, $header;
			unless ($missing{$parts[0]}) { $hash{$parts[0]} = "$parts[1]\n$_"; }
		}
	}
}
close IN1;

#print out R1
open OUT1, ">$R1File.sorted";
foreach my $entry (@master) { 
	if ($hash{$entry}) { print OUT1 "$entry $hash{$entry}\n"; }
}
close OUT1;
%hash = ();

#print missing R2
open MISS, ">$R1File.nomatch";
foreach my $entry (keys %missing) { print MISS "$entry\n"; }
close MISS;



