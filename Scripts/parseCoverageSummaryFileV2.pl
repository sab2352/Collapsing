#############################################
#
#  Parse coverage summary file
#  Use the following script to make a matrix for pca from the coverage summary files
#
# Sahar Gelfman 03/2017
#############################################

#use strict;
use warnings;
use Getopt::Long 'HelpMessage';
use Cwd;

#server
GetOptions(
  'cov_details=s' => \my $fullName,
  'sample_file=s'   => \my $fullSamplesName,
  'cluster=s'   => \my $cluster,
  'help'     =>   sub {HelpMessage(0)},
) or HelpMessage(1);

=head1 NAME

parseCoverageSummaryFileV2.pl - make a matrix for pca from the coverage details files

=head1 SYNOPSIS

parseCoverageSummaryFileV2.pl [options]

Options:

  --cov_details,-c     file path for *coverage_details.csv
  --sample_file,-s     file path for sample file
  --cluster,-u         current cluster being run
  --help,-h            Print this help

=cut

my $randomGenesFile = getcwd()."/DefaultData/randomGenesForpca.txt";
my $out=$fullName.".READYforPCA.txt";
my $log_path = getcwd()."/Data/covPCA_log/";
my $log =$log_path."parseCov_".$cluster.".log";

mkdir($log_path, 0755) unless(-d $log_path );

open (Insamples, $fullSamplesName) or die "could not find sampleFile $fullSamplesName";
open (In, $fullName) or die "could not find $fullName";
open (OUT, ">".$out) or die "could not find $out";
open (OUTLog, ">",$log) or die "could not find $log";
open (INgenes, $randomGenesFile) or die "could not find $randomGenesFile";

my %samples;
my $lineNo=0;
#Load samples and types
while (my $line = <Insamples>)
{
	my $infoLine=ProcessLine($line);	
	my (@data)= split(/\s/,$infoLine); # split by space
	$samples{$data[0]}=$data[5];
}

# loading  genes 
my %selectedGenes;
while (my $line = <INgenes>)
{
	my $infoLine=ProcessLine($line);
	$selectedGenes{$infoLine}=1;
	print OUTLog "selected gene:$infoLine\n";  
}

$lineNo=0;
#Load samples
my $header=<In>;
my $count=0;
my %sampleHash;
while (my $line = <In>)
{
	my $infoLine=ProcessLine($line);
	my (@data)= split(/,/,$infoLine); # split by space
	my $sample=$data[0];
	my $gene=$data[1];
	my $cover=$data[4];
	if (!exists($selectedGenes{$gene}))
	{
		next;
	}
	if (!exists($samples{$sample}))
	{
		next;
	}
	print OUTLog "Found gene, record no. ".$count++."\t$sample $gene \n";
	$sampleHash{$sample}{$gene}=$cover;
}

$header = "Samples";
# for each Gene (columns)
foreach my $gene (sort keys %selectedGenes)
{
	$header=$header."\t$gene"; 
}
print OUT $header."\tSample_TYPE\n";
 
#printing for each SAMPLE (lines)
foreach my $sample (sort keys %sampleHash)
{
	my $outline = "$sample";
	# for each Gene (columns)
	foreach my $gene (sort keys %selectedGenes)
	{
		if (exists($sampleHash{$sample}{$gene}))
		{
			$outline=$outline."\t".$sampleHash{$sample}{$gene}; 
		}
		else
		{
			$outline=$outline."\t0"; 
		}
	}
	# add type and print
	print OUT $outline."\t".$samples{$sample}."\n";
}
 

sub ProcessLine
{
        my $line=shift;
        $line=~s/\r\n/\n/g;
        chomp($line);
        return($line);
}
