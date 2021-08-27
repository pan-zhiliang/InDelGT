#!/usr/bin/perl -w
use strict;
use locale;
use Getopt::Long;

sub prtHelp{
        print "\nThis procedure is used for indels genotyping of F1 generation populations\n";
        print "Contact: Chunfa Tong  <tongchf\@njfu.edu.cn>\n";
        print "Usage: perl genotyping.pl [Options] <*.gt> [Options] <*.sorted.bam> [Options] <directory>\n";
        print "\n";
        print "Options:\n";
        print "	-v  <file>	the InDel genotype input file for parents or progeny\n";
        print "	-b  <file>	the sorted BAM input file of parents or progeny\n";
        print "	-o  <directory>	create a directory for storing output file of the genotyping results\n";
        print "	--help|h	help\n";

}

my $inputfile;
my $progeny_inputfile;
my $outputfile;
my $HOMOZYGOUS;
my $help;
my %hash;
my %prs;
my $samtools;
my $reference;
my ($str1,$str2);
my $number=0;
my $line_number=0;
my $length;
my (@arr,@arr1,@arr2,@arr3,@arr4);
if(@ARGV<3){
        prtHelp();
        exit;
}

GetOptions(
        "v:s"=>\$inputfile,
        "b:s"=>\$progeny_inputfile,
        "o:s"=>\$outputfile,
        "help|h" => \$help
);

if ($help){
        prtHelp();
        exit;
}

unless($inputfile){
        print STDERR "\nError: the InDel genotype input file for parents or progeny was not provided!\n\n";
        prtHelp();
        exit;
}

unless($progeny_inputfile){
        print STDERR "\nError: the sorted BAM input file of parents or progeny was not provided!\n\n";
        prtHelp();
        exit;
}

unless($outputfile){
        print STDERR "\nError: please create a directory for storing output file!\n\n";
        prtHelp();
        exit;
}

open(PRS,"../parameters.ini") || die "\nError\! This program needs a parameter file,namely \"parameters.ini\".\n\n";
while(<PRS>){
        chomp;
        if($_=~/\:/ && ($_=~/SAMTOOLS_FOLD/ or $_=~/REFERENCE_GENOME_FILE/ or $_=~/HOMOZYGOUS_DEPTH/)){
                my ($str1,$str2) = split /:/,$_;
                if($str1 =~ /\s*(\S+)\s*/){
                        $str1 = $1;
                        if($str2 =~ /\s*(\S+\s+\S+)\s*/ || $str2 =~ /\s*(\S+)\s*/){
                                $str2 = $1;
                                $prs{$str1} = $str2;
                        }
                }
        }else{
                next;
        }
}
close PRS;
if(exists($prs{"SAMTOOLS_FOLD"})){
        $samtools = $prs{"SAMTOOLS_FOLD"};
}else{
        die "Error! Please check the line of SAMTOOLS_FOLD in \"parameters.ini\".\n";
}

if(exists($prs{"REFERENCE_GENOME_FILE"})){
        $reference = $prs{"REFERENCE_GENOME_FILE"};
}else{
        die "Error! Please check the line of REFERENCE_GENOME_FILE in \"parameters.ini\".\n";
}
if(exists($prs{"HOMOZYGOUS_DEPTH"})){
        $HOMOZYGOUS = $prs{"HOMOZYGOUS_DEPTH"};
}else{
        die "Error! Please check the line of HOMOZYGOUS_DEPTH in \"parameters.ini\".\n";
}

`grep "0/0" $outputfile/$inputfile > $outputfile/$inputfile.allelloci`;
open IN, "$outputfile/$inputfile.allelloci";
while(<IN>){
	chomp;
	my @base=split/\s+/,$_;
	$length=length($base[1]);
	`$samtools tview $progeny_inputfile -p $base[0] -d T $reference > $outputfile/$inputfile-$base[0].out`;
	open FA, "$outputfile/$inputfile-$base[0].out";
	open OU, ">$outputfile/$inputfile-$base[0].trueout";
	readline FA;
	my $refline=readline FA;
	readline FA;
	my @ref=$refline=~/.{$length}/g;
	if($ref[0]=~/\*/){
		unlink "$outputfile/$inputfile-$base[0].out";
	}elsif($ref[0]!~/\*/){
		while(<FA>){
			chomp;
			if($_=~/^\s+/){
				next;
			}
			@arr = $_=~/.{$length}/g;
			if($arr[0]=~/\*/ or $arr[0]=~/[ATCG]+/){
				$number++;
			}if($arr[0]!~/\*/ && $arr[0]!~/\s+/i && $arr[0]!~/[ATCG]+/){
				$line_number++;
			}
		}
		if($number==0 && $line_number>=$HOMOZYGOUS){
			print OU "$base[0]\t$base[1]\t$base[2]\t$base[3]\n";
		}
	}
	`cat $outputfile/$inputfile-$base[0].trueout >> $outputfile/$inputfile-00.loci`;
	unlink "$outputfile/$inputfile-$base[0].trueout","$outputfile/$inputfile-$base[0].out";
	$number=0;$line_number=0
}
unlink "$outputfile/$inputfile.allelloci";
open IN, "$outputfile/$inputfile-00.loci";
while(<IN>){
        chomp;
        my @line=split/\s+/,$_;
        $hash{$line[0]}=1;
}
close IN;
open IN, "$outputfile/$inputfile";
open OU, ">$outputfile/$inputfile-1";
while(<IN>){
        chomp;
        my @line=split/\s+/,$_;
        if($_!~/0\/0/){
                print OU "$_\n";
        }elsif(exists $hash{$line[0]}){
                print OU "$_\n";
        }else{
                print OU "$line[0]\t$line[1]\t$line[2]\t_/_\n";
        }
}
undef %hash;
unlink "$outputfile/$inputfile-00.loci","$outputfile/$inputfile";
rename "$outputfile/$inputfile-1","$outputfile/$inputfile";