use strict;
use warnings;
use Digest::MD5 qw(md5_hex);

# This script will take a multiplexed Illumina FASTQ file with barcodes on the R1 (and optionally on the R2) 
# and create R1 and R2 files for each individual
# A file "summaryFile" will list their read numbers and MD5 checksums

### Hardware requirements:
# single processor
# Enough RAM for approximately two of your original files UNCOMPRESSED (i.e. a single pair of raw R1 and R2 files)
# I found that for several of our early runs, ~5Gb of RAM was needed
# Enough harddrive space for a full copy of the compressed files (split by individual) plus a full uncompressed copy of those files
# Run time is a few hours to a day, depending on the size of the data received


######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
#########  PLEASE DEFINE SOME BASIC INFORMATION THAT WE'LL NEED...

#Our data files usually have a format like this: lane8_Undetermined_L008_R1_015.fastq.tgz
#We need to split up that information so that this script can open and iterate through all files

#First, what is the filename before _R1 or _R2 and the numbers? (Usually the same, sometimes not...)
my $baseNameR1 = '3383_Heather';
my $baseNameR2 = '3383_Heather';

#How many R1/R2 file sets are there? You'll need to look inside the run folder and see what they count up to.
my $setNum = '1';

#What is the compressed file extension for this run? indicate 'tar.gz' or 'gz' or 'tgz' or 'NONE'  (any others will require altering this script)
my $ext = 'NONE';

#What is the full location name of the run folder that will be processed? Starting with /home (no trailing /)
my $run = '/home/kdlugosch/data_temp/Heather_run_02Dec2016';

#What is the full location name of a barcode file that lists by row:
#the Sample ID, its P1 barcode, and its P2 revcomp barcode? (three tab delimited columns, no header line)
#NOTE!!! The P2 (R2) barcode must be the reverse complement to the sequence listed on the spreadsheet/tube!!!!!!
#For samples with no barcode on the R2, use the word NONE (all caps) in place of a code.
#Do not make this file in a word processor - use a simple text editor that won't put in hidden characters.
my $barcodeKey = '/home/kdlugosch/data_temp/Heather_run_02Dec2016/barcodes';

#What is the full location where you want the resulting folder of UNCOMPRESSED FASTQ files for each individual?
my $outLocation = '/home/kdlugosch/data_temp/Heather_run_02Dec2016';
#A new folder will be made in that location. What name do you want for this?
my $outDir = "$outLocation/testout2";

#Set the allowed number of mismatches between expected and observed barcodes
my $allowed_mismatches = 1;

#run: perl Illumina_splitting_script_[version].pl

######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 





######### NOW TO GET THINGS READY
#Check that the filenames that we have entered are correct (-e askes if something exists)
unless ($ext eq 'NONE') {
	unless (-e "$run/$baseNameR1\_R1_001.fastq.$ext") { die "\n\nError: file $run/$baseNameR1\_R1_001.fastq.$ext does not exist!\n\n"; }
	unless (-e "$run/$baseNameR2\_R2_001.fastq.$ext") { die "\n\nError: file $run/$baseNameR2\_R2_001.fastq.$ext does not exist!\n\n"; }
}
else {
	unless (-e "$run/$baseNameR1\_R1_001.fastq") { die "\n\nError: file $run/$baseNameR1\_R1_001.fastq does not exist!\n\n"; }
	unless (-e "$run/$baseNameR2\_R2_001.fastq") { die "\n\nError: file $run/$baseNameR2\_R2_001.fastq does not exist!\n\n"; }
}


unless (-e "$barcodeKey") { die "\n\nError: file $barcodeKey does not exist!\n\n"; }

#So that we can work with the barcode information, read in the list of barcodes
my (%barcodes, $codeLen);
my (%R1list, %R2list);
my $noR2code = 'N'; #initialize a record of whether any R2 are listed as 'NONE'
open CODES, "<$barcodeKey";
foreach (<CODES>) {
	chomp $_;
	my @cols = split /\t/, $_; #cols should be [0]ID [1]R1 barcode [2] R2 barcode
	if (($cols[0]) && ($cols[1]) && ($cols[2])) {
		$codeLen = length($cols[1]);
		$cols[1] = uc $cols[1]; $cols[2] = uc $cols[2];
		$R1list{$cols[1]}++; $R2list{$cols[2]}++; 
		$barcodes{"$cols[1].$cols[2]"} = $cols[0];
		if ($cols[2] eq 'NONE') { $noR2code = 'Y'; }
	}
	else { die "\n\nError: Barcode file does not appear to have three tab deliminated columns.\n\n"; }
}
close CODES;

#Make new folders for the run analysis
unless (-e "$outDir") { system ("mkdir $outDir"); }
#system ("cp splitBC_Sno_noAmbig.pl $outDir/");

#Make the barcode files for the splitting program
#open R1CODE, ">$outDir/R1_codes";
#foreach my $code (keys %R1list) { print R1CODE "$code\t$code\n"; }
#close R1CODE;

my $R2count = 0; my $R2split = 'Y';
#open R2CODE, ">$outDir/R2_codes";
foreach my $code (keys %R2list) { unless ($code eq 'NONE') { $R2count++; }}
#print R2CODE "$code\t$code\n";
#close R2CODE;
if ($R2count == 0) { $R2split = 'N'; } #If all of the R2 are NONE, then there is no splitting of that file
	
#prep output files for sequences that don't match a barcode set
open NOR1, ">$outDir/noMatch_R1.fastq"; close NOR1;
open NOR2, ">$outDir/noMatch_R2.fastq"; close NOR2;
	
	


######### NOW LET'S PROCESS OUR DATA
my (%readCount, %R1count, %R2count);
#For each set of files, move them to the new folder and split them
my ($R1INFILE, $R2INFILE);
for (my $j = 1; $j <= $setNum; $j++) {

	print "\nNow processing file pair $j for $run";
	
	############ DEFINE THE UNPACKED FILENAMES
	if ($j < 10) { 
		$R1INFILE = "$baseNameR1\_R1_00$j.fastq"; #single digit file count needs an extra zero in the name
		$R2INFILE = "$baseNameR2\_R2_00$j.fastq"; #single digit file count needs an extra zero in the name
	}
	else {
		$R1INFILE = "$baseNameR1\_R1_0$j.fastq"; #double digit file count needs one less zero in the name
		$R2INFILE = "$baseNameR2\_R2_0$j.fastq"; #double digit file count needs one less zero in the name
	}
	
	
	############ R1 barcode reading
	#Unpack the R1 into the working folder
	my $R1openFile;
	if (($ext eq 'tgz') || ($ext eq 'tar.gz')) { system ("cd $outDir/; tar -xzf $run/$R1INFILE.$ext"); $R1openFile = "$outDir/$R1INFILE";}
	elsif ($ext eq 'gz') { system ("gunzip -c $run/$R1INFILE.$ext >$outDir/$R1INFILE"); $R1openFile = "$outDir/$R1INFILE";}
	elsif ($ext eq 'NONE') { $R1openFile = "$run/$R1INFILE"; }#system ("cp $run/$R1INFILE $outDir/$R1INFILE"); }  #TODO remove this?
	else { die "\n\nError: Compressed file type '$ext' is not a recognized option for this script\n\n"; }
	
	#Read the unpacked R1
	print "\nReading R1-$j...";
	my (%R1hash, $header, $headerAll, $seq, $qual) = ();
	my $i = 0; #start an index to count lines - four lines to an entry
	open R1, "<$R1openFile";
	while (<R1>) {
		$i++; #keep track of lines, counting to 4 and then starting over
		if ($i == 1) { $headerAll = $_; my @parts = split /\s/, $headerAll; $header = $parts[0];} #header line
		elsif ($i == 2) { $seq = $_; } #DNA line
		elsif ($i == 3) { chomp $_; unless ($_ eq '+') { die "\n\nError: FASTQ format for $outDir/$R1INFILE is not reading correctly. Where '+' expected, there is $_\n\n"; }} #plus line
		elsif ($i == 4) { #Qual line, end of entry
			my $codeArea = substr($seq, 0, $codeLen);
			$R1hash{$header} = "$codeArea"; #record the header and barcode
			($headerAll, $header, $seq, $qual) = (); #reset everything
			$i = 0
		} 
	}
	close R1;
	

	############ R2 reading and splitting
	print "\nReading and splitting R2-$j...\n";
	
	#Unpack the R2 into the working folder
	my $R2openFile;
	if (($ext eq 'tgz') || ($ext eq 'tar.gz')) { system ("cd $outDir/; tar -xzf $run/$R2INFILE.$ext"); $R2openFile = "$outDir/$R2INFILE";}
	elsif ($ext eq 'gz') { system ("gunzip -c $run/$R2INFILE.$ext >$outDir/$R2INFILE"); $R2openFile = "$outDir/$R2INFILE"; }
	elsif ($ext eq 'NONE') { $R2openFile = "$run/$R2INFILE"; }
	else { die "\n\nError: Compressed file type '$ext' is not a recognized option for this script\n\n"; }
	
	#read in R2 file
	my ($id);
	($headerAll, $header, $seq, $qual) = ();
	$i = 0; #start an index to count lines - four lines to an entry

	open R2, "<$R2openFile";
	while (<R2>) { #read in sequences
		$i++; #keep track of lines, counting to 4 and then starting over
		#chomp $_; 
		if ($i == 1) { $headerAll = $_; my @parts = split /\s/, $headerAll; $header = $parts[0];} #header line
		elsif ($i == 2) { $seq = $_; } #DNA line
		elsif ($i == 3) { chomp $_; unless ($_ eq '+') { die "\n\nError: FASTQ format for $outDir/$R2INFILE is not reading correctly. Where '+' expected, there is $_\n\n"; }} #plus line
		elsif ($i == 4) { #Qual line
			$qual = $_;
			unless ($R1hash{$header}) { die "\n\nError: Header '$header' in the R2 file doesn't have a match in the R1 list\n\n";}
			my $codeArea = substr($seq, 0, $codeLen);

			#figure out the R1 barcode match if any 
			my $R1codeAssign;
			if ($R1list{$R1hash{$header}}) { $R1codeAssign = "$R1hash{$header}"; } #there is an exact match to an expected barcode
			else { #there is no exact match, so see if there is one within in 1 bp mismatch
				$R1codeAssign = match_sequences($R1hash{$header},\%R1list);
				#if this is assigned to be NONE, then the sequence is entirely unknown
				if ($R1codeAssign eq 'NONE'){ $id = 'noMatch'; }
			}
			$R1count{$R1codeAssign}++;

			#figure out the R2 barcode match if any  
			my $R2codeAssign;
			if ($R1codeAssign eq 'NONE'){ $R2codeAssign = 'NONE'; }
			else {
			if ($R2split eq 'N') { $R2codeAssign = 'NONE'; } #there are no R2 barcodes in this dataset
			elsif ($R2list{$codeArea}) { $R2codeAssign = "$codeArea"; } #there is an exact match to an expected barcode
			else { #there is no exact match, so see if there is one within in 1 bp mismatch
				$R2codeAssign = match_sequences($codeArea,\%R2list);
			}}
			$R2count{$R2codeAssign}++;

			#assign the individual
			unless ($R1codeAssign eq 'NONE'){
				my $codeCombo = "$R1codeAssign.$R2codeAssign";
				if ($barcodes{$codeCombo}) { $id = $barcodes{$codeCombo}; }
				else { $id = 'noMatch'; }
			}

			#print out R2 data to appropriate files
			my $R2printArea; 
			if ($id eq 'noMatch') { $R2printArea = $seq; } 
			#trim the barcode from both the sequence and qual lines
			else { $R2printArea = substr($seq, $codeLen); $qual = substr($qual, $codeLen);}
			open OUT2, ">>$outDir/$id\_R2.fastq"; print OUT2 "$headerAll$R2printArea+\n$qual"; close OUT2;
	
			#roundup
			$readCount{$id}++;
			$R1hash{$header} = $id; #replace barcode entry for R1 with its individual ID
			($headerAll,$header, $seq, $qual, $id) = (); #reset everything
			$i = 0
		} 
	}
	close R2;
	unless ($ext eq 'NONE') { system ("rm $outDir/$R2INFILE"); }

	############ R1 splitting
	print "\nSplitting R1-$j...\n";
	($header, $headerAll, $seq, $qual) = ();
	$i = 0; #start an index to count lines - four lines to an entry
	my $destination;
	open R1, "<$R1openFile";
	while (<R1>) {
		$i++; #keep track of lines, counting to 4 and then starting over
		if ($i == 1) { 
			my @parts = split /\s/, $_; $header = $parts[0]; #header line
			$destination = $R1hash{$header};
			open OUT1, ">>$outDir/$destination\_R1.fastq"; print OUT1 "$_";
		}
		elsif ($i == 2) { #DNA line
			$seq = $_;
			my $R1printArea; if ($destination eq 'noMatch') { $R1printArea = $seq;} else {$R1printArea = substr($seq, $codeLen);}
			print OUT1 "$R1printArea"; 
		} 
		elsif ($i == 3) { #plus line
			if ($_ =~ m/^\+/) { print OUT1 "$_"; }
			else { die "\n\nError: FASTQ format for $outDir/$R1INFILE is not reading correctly. Where '+' expected, there is $_\n\n"; }
		} 
		elsif ($i == 4) { #Qual line, end of entry
			my $qualPrint; if ($destination eq 'noMatch') { $qualPrint = $_;} else {$qualPrint = substr($_, $codeLen);}
			print OUT1 "$qualPrint";
			close OUT1;
			($headerAll, $destination) = (); #reset everything
			$i = 0
		} 
	}
	close R1;
	%R1hash = ();
	unless ($ext eq 'NONE') { system ("rm $outDir/$R1INFILE"); }
}


#Report stats on barcode matches in R1count, R2count
print "\nReporting barcode matches in split_log...\n";
open COUNTS, ">$outDir/split_log";
print COUNTS "R1 Barcode\tcount\n";
foreach my $entry (keys %R1count){ print COUNTS "$entry\t$R1count{$entry}\n"; }
print COUNTS "\nR2 Barcode\tcount\n";
foreach my $entry (keys %R2count){ print COUNTS "$entry\t$R2count{$entry}\n"; }
close COUNTS;


#Now that files have been split into a set for each individual by the subroutine, we want to record their final checksum
print "\nCalculating checksums...\n";
open HASH, ">$outDir/hashFile.md5";
open INFO, ">$outDir/summaryFile";
print INFO "Individual\tReadNumber\tR1checksum\tR2checksum\n";
my ($R1check, $R2check);
foreach my $sample (keys %readCount){
  unless ($sample eq 'noMatch'){
	open R1FINAL, "<$outDir/$sample\_R1.fastq"; my $R1check = md5_hex(<R1FINAL>); close R1FINAL;
	open R2FINAL, "<$outDir/$sample\_R2.fastq"; my $R2check = md5_hex(<R2FINAL>); close R2FINAL;
	print INFO "$sample\t$readCount{$sample}\t$R1check\t$R2check\n";
	print HASH "$R1check  $sample\_R1.fastq\n$R2check  $sample\_R2.fastq\n";
  }
  else {print INFO "$sample\t$readCount{$sample}\tNA\tNA\n";} #only readcount information for noMatch file
}
close INFO;
close HASH;

print "\nProcessing of $run is complete!\n\n\n";



#modified from SplitBC_Sno_noAmbig.pl
sub match_sequences {

	my $sequence_fragment = shift; #This is the seq barcode fragment, the first element fed in
	my $barcodeHashRef = shift; #This the R1 or R2 list ref, the second element fed in

	# split accotding to barcodes
	my $best_barcode_mismatches_count = $codeLen; 
	my $best_barcode_ident = undef;

	#Try all barcodes, find the one with the lowest mismatch count
	my %misCount; #from KMD = let's keep track of how many times a mismatch count is seen (in case there are two tied best matches)
	foreach my $barcoderef (keys %{$barcodeHashRef}) { 
		my $ident = $barcoderef; my $barcode = $barcoderef;
		my $mm = mismatch_count($sequence_fragment, $barcode); 

		#KMD Note: The following single if loop is too simple. It simply records the lowest mismatch score seen 
		#and assoc barcode identity. If there are two equivalent matches, this is not recorded anywhere and whichever 
		#was seen last will be the assignment. 
		# First, we need a way to record how many times a mismatch count is seen. A count hash should do the trick.
		$misCount{$mm}++; #New from KMD - records a count for every mismatch value seen
			
		if ( $mm < $best_barcode_mismatches_count ) {
			$best_barcode_mismatches_count = $mm ;
			$best_barcode_ident = $ident ;
		}
	}
		
	if ($misCount{$best_barcode_mismatches_count}>1) {$best_barcode_mismatches_count = $codeLen;} #New from KMD - this resets 
	#the mismatch count to be too high if there is an ambiguous match, i.e. more than one match with the lowest mismatch count.
		
	if ( (!defined $best_barcode_ident) || ($best_barcode_mismatches_count>$allowed_mismatches)){
		$best_barcode_ident = 'NONE'; 
	}

	return ($best_barcode_ident);
}

#Quickly calculate hamming distance between two strings
#
#NOTE: Strings must be same length.
#      returns number of different characters.
#see  http://www.perlmonks.org/?node_id=500235
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }


