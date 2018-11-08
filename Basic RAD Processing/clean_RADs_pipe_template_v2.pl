#!usr/bin/perl
use strict;
use warnings;
use Cwd;

# This script will take one or more folders of split RAD files (fastq format) and put them through SnoWhite cleaning
# It assumes:
#	1) Inside the given folder(s), all files of interest start with CESO (XXX change on line 4)
#	2) This script is inside the snowhite run folder (in the same location as the snowhite_2.0.3.pl script, and snowhite is ready to go
#	3) See line 20 for setting up source and destination folders
#	4) See lines 59 and 87 XXX for snowhite commands and update as needed
#
# run:	perl clean_RADs_pipe.pl
#
##############################


#Set up the source and destination folders
my %folders;
#FOR YOU TODO Enter the full source and destination folder locations in the following format
#$folders{'locationOfInputFiles'} = 'locationOfOutputfiles'; NOTE that the output folder will be created if it doesn't already exist
#here is an example of a real one - USE # to comment this out of the code or use the line for your own information!
$folders{'/data/2012_11_15_Illumina_100bpPE_split_exampleFolder'} = '/data/YST_RADs_cleanedTrimmed/2012_11_15_Illumina_100bpPE'; #run1

#Set up the trimming lengths
#recognition site left on split reads: PstI:5bp; EcoRI:5bp; MseI:3bp
my $R1trim = 5; #length of R1 (forward) enzyme recognition site to trim
my $R2trim = 3; #length of R2 (reverse) enzyme recognition site to trim
my $R1length = 76; #length of R1 (forward) fragment to keep [76bp was used in YST MEC paper, a 13bp loss]
my $R2length = 70; #length of R2 (reverse) fragment to keep [these tend to loss quality more quickly than R1, this is a guess] #TODO


############SCRIPT STARTS HERE
my $dir = cwd();
my $filename;
#go though each folder
foreach my $sourceDir (keys %folders) {
	#harvest the files; they have form 'CESO_Gr18-24-5_R1.fastq'
	#in this case there are some extraneous files that have .fastq endings, so we want only those that start with CESO
	my @files = `ls $sourceDir/CESO_*`; #XXX make sure the correct name is here
	#my @files = `ls $sourceDir/*`;

	#foreach of those files
	foreach my $file (@files) {
		#remove the path
		chomp $file;
		my @pathParts = split /\//, $file;
		$filename = $pathParts[-1]; 

		#clean up from the last run
		system ("rm -rf temp");

		#start processing this file
		if ($filename =~ m/R1\.fastq$/) {
			system ("cp $file $dir");
			
			#run snowhite
			print "\n\nCleaning file <$file>";
			system ("perl snowhite_2.0.3.pl -f $filename -v RAD_primers_Dec2011 -s RAD_primers_Dec2011  -g T -Q 10 -D T -L T -m $R1length -E $R1trim -c 5 -p 3 -o temp"); 
			# -g T	deletes temp files
			# -Q 10	min phred score, 10 is recommended
			# -D T	run TagDust
			# -L T	run SeqClean
			# -m 76	min sequence length
			# -E 5	5' end clipping
			# -c 5	length of end clipping
			# -p 3	processor number

			
			#trim the sequences to length (snowhite minimum should mean little read loss here)
			print "=====\n\nTrimming cleaned file <temp/FinalOutput/temp.clean>"; &trim($R1length);

			#copy the file to destination
			#if the destination folder doesn't exist already, make it
			unless (-e "$folders{$sourceDir}") { system("mkdir $folders{$sourceDir}"); }
			system ("mv $filename.clean.trim_$R1length $folders{$sourceDir}/");
			
			system ("rm $filename");
			print "\nDone with file <$file>\n\n";

		}
		elsif ($filename =~ m/R2\.fastq$/) {
			system ("cp $file $dir");
			
			#run snowhite
			print "\n\nCleaning file <$file>";
			system ("perl snowhite_2.0.3.pl -f $filename -v RAD_primers_Dec2011 -s RAD_primers_Dec2011  -g T -Q 10 -D T -L T -m $R2length -E $R2trim -c 3 -p 3 -o temp"); 

			#trim the sequences to length
			print "=====\n\nTrimming cleaned file <temp/FinalOutput/temp.clean>"; &trim($R2length);

			#copy the file to destination
			#if the destination folder doesn't exist already, make it
			unless (-e "$folders{$sourceDir}") { system("mkdir $folders{$sourceDir}"); }
			system ("mv $filename.clean.trim_$R2length $folders{$sourceDir}/");
			
			system ("rm $filename");
			print "\nDone with file <$file>\n\n";
		}
		
	}	
}
print "\nDone will all files!\n Note that the snowhite output from the last file has been left for inspection in the 'temp' folder.\n\n";


sub trim {
	my $length = shift;

	#initiate output file
	open OUT, ">$filename.clean.trim_$length"; 

	##Read in file
	my $header = 'Expecting FASTA file';
	open IN, "<temp/FinalOutput/temp.clean" || die "\n\tError: Can't open file <temp/FinalOutput/temp.clean>\n\n";
	my $sequence = ();
	foreach(<IN>) { 
		chomp $_; 
		unless ($_) { next; }  # skip empty lines
   		if ($_ =~ m/^>/) { #line is a header
			#if a sequence has been accumulated from last entry, process the last entry
			if ($sequence) {
				if (length($sequence) >= $length) {
					my $trimmed = substr($sequence,0,$length); #trim 
					print OUT "$header\n$trimmed\n";
				}
			}
			#store the new header and clear the sequence
			$header = $_;
			$sequence = ();
		}
		else { $sequence .= $_; }

	}
	close IN;
	close OUT;

}

