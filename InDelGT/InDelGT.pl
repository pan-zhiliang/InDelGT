#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Parallel::ForkManager;

sub prtHelp{
	print "\nInDelGT is an integrated pipeline for extracting InDel genotypes in a hybrid population with next generation DNA sequencing data.\n";
	print "Contact: Chunfa Tong <tongchf\@njfu.edu.cn>\n";
	print "Version: 1.0\n";
	print "Usage: perl InDelGT.pl [option] <type> -o <directory>\n";
	print "\n";
	print "Options:\n";
	print "	-p  <type>      the type of population: CP;BC1;BC2;F2 (default: CP)\n";
	print "	-o <directory>	create a directory for storing output file of the InDel genotyping results\n";
	print "	--help|h	help\n";
}

my %prs;
my ($str1,$str2);
my $outDirName;
my $population='CP';
my $help;
my $np=0;
if(@ARGV<1){
	prtHelp();
	exit;
}

GetOptions(
	"p:s"=>\$population,
	"o:s"=>\$outDirName,
	"help|h"=>\$help
);


if($help){
	prtHelp();
	exit;
}

unless($population eq 'CP' or $population eq 'BC1' or $population eq 'BC2' or $population eq 'F2'){
        print STDERR "\nError: option '-p' only be set to 'CP' or 'BC1' or 'BC2' or 'F2'!\n\n";
        prtHelp();
        exit;
}

unless($outDirName){
	print STDERR "\nError: please create a directory for storing output file!\n\n";
	prtHelp();
	exit;
}

open(PRS,"<parameters.ini") || die "\nError\! This program needs a parameter file,namely \"parameters.ini\".\n\n";
while(<PRS>){
	chomp;
	if(/:/){
		($str1,$str2) = split /:/,$_;
		if($str1 =~ /\s*(\S+)\s*/){
			$str1 = $1;
			if($str2 =~ /\s*(\S+\s+\S+)\s*/ || $str2 =~ /\s*(\S+)\s*/){
				$str2 = $1;
				$prs{$str1} = $str2;
				if($str1=~/^\s*PROGENY\d+/){
					$np++;
				}
			}
		}
	}
}
close PRS;

my ($bwa,$samtools,$bcftools,$INDELGT,$datafold,$reference,$parent1,$parent2,$progeny_number,$MAPQ,$threads);

if(exists($prs{"BWA_FOLD"})){
	$bwa = $prs{"BWA_FOLD"};
}else{
	die "Error! Please check the line of BWA_FOLD in \"parameters.ini\".\n";
}

if(exists($prs{"SAMTOOLS_FOLD"})){
	$samtools = $prs{"SAMTOOLS_FOLD"};
}else{
	die "Error! Please check the line of SAMTOOLS_FOLD in \"parameters.ini\".\n";
}

if(exists($prs{"BCFTOOLS_FOLD"})){
	$bcftools = $prs{"BCFTOOLS_FOLD"};
}else{
	die "Error! Please check the line of BCFTOOLS_FOLD in \"parameters.ini\".\n";
}

if(exists($prs{"DATAFOLD"})){
	$datafold = $prs{"DATAFOLD"};
}else{
	die "Error! Please check the line of DATAFOLD in \"parameters.ini\".\n";
}

if(exists($prs{"INDELGT_FOLD"})){
	$INDELGT = $prs{"INDELGT_FOLD"};
}else{
	die "Error! Please check the line of INDELGT_FOLD in \"parameters.ini\".\n";
}

if(exists($prs{"REFERENCE_FILE"})){
	$reference = $prs{"REFERENCE_FILE"};
}else{
	die "Error! Please check the line of REFERENCE_FILE in \"parameters.ini\".\n";
}

if(exists($prs{"PARENT1"})){
        $parent1 = $prs{"PARENT1"};
}else{
        die "Error! Please check the line of PARENT1 in \"parameters.ini\".\n";
}

if(exists($prs{"PARENT2"})){
        $parent2 = $prs{"PARENT2"};
}else{
        die "Error! Please check the line of PARENT2 in \"parameters.ini\".\n";
}

if(exists($prs{"PROGENY_NUMBER"})){
        $progeny_number = $prs{"PROGENY_NUMBER"};
}else{
        die "Error! Please check the line of PROGENY_NUMBER in \"parameters.ini\".\n";
}

if(exists($prs{"MAPQ"})){
        $MAPQ = $prs{"MAPQ"};
}else{
        die "Error! Please check the line of MAPQ in \"parameters.ini\".\n";
}

if(exists($prs{"THREADS"})){
        $threads = $prs{"THREADS"};
}else{
        die "Error! Please check the line of THREADS in \"parameters.ini\".\n";
}

####################################Identification of male parent indels###############################################

my $ref="$datafold/$reference";
my @parent1_fq=split/\s+/,$parent1;
my @parent2_fq=split/\s+/,$parent2;
my @refGenome=split/\//,$ref;
my $parent1_fq_file1="$datafold/$parent1_fq[0]";
my $parent1_fq_file2="$datafold/$parent1_fq[1]";
my $parent2_fq_file1="$datafold/$parent2_fq[0]";
my $parent2_fq_file2="$datafold/$parent2_fq[1]";
if(-d $outDirName){
        `rm -r $outDirName`;
}if(!-d $outDirName){
	mkdir $outDirName, or die "Error: can't create directory '$outDirName' : $!" unless( -d $outDirName);
	chdir $outDirName or die "Error: can't cd to directory '$outDirName' : $!";
}
`cp $ref $INDELGT/$outDirName`;
if(!-e $ref){
	print STDERR "Error: The reference genome was not prepared. Please prepare the reference genome!\n\n";
	exit;
}if(-e $reference){
	if($bwa=~/bwa-mem2/){
		`$bwa/bwa-mem2 index $reference`;
	}else{
		`$bwa/bwa index $reference`;
	}
}
if((!-e $parent1_fq_file1) or (!-e $parent1_fq_file2)){
	print STDERR "Error: The pair-end reads of the parent1 were not prepared. Please prepare the pair-end reads of the parent1!\n\n";
	exit;
}if((-e $parent1_fq_file1) or (-e $parent1_fq_file2)){
	if($bwa=~/bwa-mem2/){
		`$bwa/bwa-mem2 mem -t $threads $reference $parent1_fq_file1 $parent1_fq_file2 > parent1.sam`;
	}else{
		`$bwa/bwa mem -t $threads $reference $parent1_fq_file1 $parent1_fq_file2 > parent1.sam`;
	}
}
`$samtools/samtools sort -@ $threads -n -o parent1.sorted.sam parent1.sam`;
if(-e "parent1.sorted.sam"){
	unlink ("parent1.sam");
}
`$samtools/samtools fixmate -@ $threads -m parent1.sorted.sam parent1.sorted.fixmate`;
if(-e "parent1.sorted.fixmate"){
	unlink ("parent1.sorted.sam");
}
`$samtools/samtools sort  -@ $threads parent1.sorted.fixmate -o parent1.sorted.fixmate1`;
if(-e "parent1.sorted.fixmate1"){
	unlink ("parent1.sorted.fixmate");
}
`$samtools/samtools markdup -@ $threads  -r parent1.sorted.fixmate1 parent1.sorted.markdup`;
if(-e "parent1.sorted.markdup"){
	unlink ("parent1.sorted.fixmate1");
}
`$samtools/samtools view -@ $threads -b -S parent1.sorted.markdup -q $MAPQ > parent1.bam`;
if(-e "parent1.bam"){
	unlink ("parent1.sorted.markdup");
}
`$samtools/samtools sort -@ $threads parent1.bam > parent1.sorted.bam`;
if(-e "parent1.sorted.bam"){
	unlink ("parent1.bam");
}
`$samtools/samtools index parent1.sorted.bam`;
`$bcftools/bcftools mpileup -Obuzv -a  AD,INFO/AD -f $reference parent1.sorted.bam > parent1.bcf`;
`$bcftools/bcftools call -m -v -f gq parent1.bcf > parent1.vcf`;

####################################Identification of female parent indels###############################################

if((!-e $parent2_fq_file1) or (!-e $parent2_fq_file2)){
	print STDERR "Error: The pair-end reads of the parent2 were not prepared. Please prepare the pair-end reads of the parent2!\n\n";
	exit;
}if((-e $parent2_fq_file1) or (-e $parent2_fq_file2)){
	if($bwa=~/bwa-mem2/){
		`$bwa/bwa-mem2 mem -t $threads $reference $parent2_fq_file1 $parent2_fq_file2 > parent2.sam`;
	}else{
		`$bwa/bwa mem -t $threads $reference $parent2_fq_file1 $parent2_fq_file2 > parent2.sam`;
	}
}
`$samtools/samtools sort -@ $threads -n -o parent2.sorted.sam parent2.sam`;
if(-e "parent2.sorted.sam"){
	unlink ("parent2.sam");
}
`$samtools/samtools fixmate -@ $threads -m parent2.sorted.sam parent2.sorted.fixmate`;
if(-e "parent2.sorted.fixmate"){
	unlink ("parent2.sorted.sam");
}
`$samtools/samtools sort  -@ $threads parent2.sorted.fixmate -o parent2.sorted.fixmate1`;
if(-e "parent2.sorted.fixmate1"){
	unlink ("parent2.sorted.fixmate");
}
`$samtools/samtools markdup -@ $threads  -r parent2.sorted.fixmate1 parent2.sorted.markdup`;
if(-e "parent2.sorted.markdup"){
	unlink ("parent2.sorted.fixmate1");
}
`$samtools/samtools view -@ $threads -b -S parent2.sorted.markdup -q $MAPQ > parent2.bam`;
if(-e "parent2.bam"){
	unlink ("parent2.sorted.markdup");
}
`$samtools/samtools sort -@ $threads parent2.bam > parent2.sorted.bam`;
if(-e "parent2.sorted.bam"){
        unlink ("parent2.bam");
}
`$samtools/samtools index parent2.sorted.bam`;
`$bcftools/bcftools mpileup -Obuzv -a  AD,INFO/AD -f $reference parent2.sorted.bam > parent2.bcf`;
`$bcftools/bcftools call -m -v -f gq parent2.bcf > parent2.vcf`;

####################################The end of identification of parents indels#############################################

##########################################Indels genotyping of the parents##################################################
chdir "$INDELGT/pls";
`perl parent_genotyping.pl -v parent1.vcf -b parent2.bcf -o $outDirName`;
`perl parent_genotyping.pl -v parent2.vcf -b parent1.bcf -o $outDirName`;
`perl homozygote_filter.pl -v parent1_parent2.parentgt -b parent1.sorted.bam -o $outDirName`;
`perl homozygote_filter.pl -v parent2_parent1.parentgt -b parent2.sorted.bam -o $outDirName`;
chdir "$INDELGT/$outDirName";
####################################The end of indels genotyping of the parents#############################################

##########################################Indels genotyping of the progeny##################################################

my $i;
my $str0;
my $flag;
my @prg1fq;
my @prg2fq;
my ($fqprfx0,@fqprfx);
for($i = 1;$i <= $np;$i++){
	$str0 = "PROGENY$i";
	if(exists($prs{$str0})){
		($prg1fq[$i-1],$prg2fq[$i-1]) = split /\s+/,$prs{$str0};
		$str1 = "$datafold\/$prg1fq[$i-1]";
		$str2 = "$datafold\/$prg2fq[$i-1]";
		unless(-e $str1){
			die "Error\! The first fastq file of progeny $i is wrong or does not exist.\n";
		}
		unless(-e $str2){
			die "Error\! The second fastq file of progeny $i is wrong or does not exist.\n";
		}

		($flag,$fqprfx0) = &prefixfq0($prg1fq[$i-1],$prg2fq[$i-1]);

		if($flag == 1){
			push @fqprfx,$fqprfx0;
		}else{
			die "Error: The file names of $prg1fq[$i-1] and $prg2fq[$i-1] are not consistent!\n";
			exit;
		}

	}else{
		die "Error! Please check the line of PROGENY$i in \"parameters.ini\".\n";
	}
	open IN, "<$INDELGT/parameters.ini";
	while(<IN>){
		chomp;
		if($_=~/PROGENY$i\:/){
			my @id=split/\:/,$_;
			if($id[1]=~/\s*(\S+\s+\S+)\s*/){
				my @fq=split/\s+/,$1;
				if($bwa=~/bwa-mem2/){
					`$bwa/bwa-mem2 mem -t $threads $reference $datafold/$fq[0] $datafold/$fq[1] > $fqprfx[$i-1]_progeny.sam`;
				}else{
					`$bwa/bwa mem -t $threads $reference $datafold/$fq[0] $datafold/$fq[1] > $fqprfx[$i-1]_progeny.sam`;
				}
				`$samtools/samtools sort -@ $threads -n -o $fqprfx[$i-1]_progeny.sorted.sam $fqprfx[$i-1]_progeny.sam`;
				if(-e "$fqprfx[$i-1]_progeny.sorted.sam"){
					unlink ("$fqprfx[$i-1]_progeny.sam");
				}
				`$samtools/samtools fixmate -@ $threads -m $fqprfx[$i-1]_progeny.sorted.sam $fqprfx[$i-1]_progeny.sorted.fixmate`;
				if(-e "$fqprfx[$i-1]_progeny.sorted.fixmate"){
					unlink ("$fqprfx[$i-1]_progeny.sorted.sam");
				}
				`$samtools/samtools sort  -@ $threads -o $fqprfx[$i-1]_progeny.sorted.fixmate1 $fqprfx[$i-1]_progeny.sorted.fixmate`;
				if(-e "$fqprfx[$i-1]_progeny.sorted.fixmate1"){
					unlink ("$fqprfx[$i-1]_progeny.sorted.fixmate");
				}
				`$samtools/samtools markdup -@ $threads  -r $fqprfx[$i-1]_progeny.sorted.fixmate1 $fqprfx[$i-1]_progeny.sorted.markdup`;
				if(-e "$fqprfx[$i-1]_progeny.sorted.markdup"){
					unlink ("$fqprfx[$i-1]_progeny.sorted.fixmate1");
				}
				`$samtools/samtools view -@ $threads -b -S $fqprfx[$i-1]_progeny.sorted.markdup -q $MAPQ > $fqprfx[$i-1]_progeny.bam`;
				if(-e "$fqprfx[$i-1]_progeny.bam"){
				unlink ("$fqprfx[$i-1]_progeny.sorted.markdup");
				}
				`$samtools/samtools sort -@ $threads $fqprfx[$i-1]_progeny.bam > $fqprfx[$i-1]_progeny.sorted.bam`;
				if(-e "$fqprfx[$i-1]_progeny.sorted.bam"){
					unlink ("$fqprfx[$i-1]_progeny.bam");
				}
				`$samtools/samtools index $fqprfx[$i-1]_progeny.sorted.bam`;
			}
		}else{
			next;
		}
	}
	close IN;
}

my $pm = new Parallel::ForkManager($threads);
foreach my $Sorted_Name (glob("*_progeny.sorted.bam")){
	$Sorted_Name=~s/\.sorted.bam//g;
	my $pid=$pm-> start and next;
	`$bcftools/bcftools mpileup -Obuzv -a  AD,INFO/AD -f $reference $Sorted_Name\.sorted.bam > $Sorted_Name.bcf`;
	$pm->finish;
}
$pm->wait_all_children;

my @Sorted_Name1=glob("*_progeny.sorted.bam");
chdir "$INDELGT/pls";
foreach my $bam (@Sorted_Name1) {
	$bam=~s/\.sorted.bam//g;
	my $pid=$pm-> start and next;
	`perl progeny_genotyping.pl -v parent1.vcf -b $bam\.bcf -o $outDirName`;
	`perl progeny_genotyping.pl -v parent2.vcf -b $bam\.bcf -o $outDirName`;
	`perl homozygote_filter.pl -v $bam\_parent1.gt -b $bam\.sorted.bam -o $outDirName`;
	`perl homozygote_filter.pl -v $bam\_parent2.gt -b $bam\.sorted.bam -o $outDirName`;
	$pm->finish;
}
$pm->wait_all_children;
chdir "$INDELGT/$outDirName";
unlink glob"$reference*";
chdir "$INDELGT/pls";
unless($population eq "BC1" or $population eq "BC2"){
	`perl segregation_ratio.pl -v parent2.vcf -p $population -d $outDirName -i $progeny_number`;
}
chdir "$INDELGT/$outDirName";
if(-e "parent2_abxaa.txt"){
	`mv parent2_abxaa.txt aaxab.txt`;
}
if(-e "parent2_abxcc.txt"){
	open IN, "parent2_abxcc.txt";
	open OU, ">aaxbc.txt";
	while(<IN>){
		chomp;
		if($_=~/Position/){
			$_=~s/ac/ab/g;
			$_=~s/bc/ac/g;
			print OU "$_\n";
		}else{
			my @line=split/\s+/,$_;
			my $aabcid=shift@line;
			my $aabcgt1=shift@line;
			my @aabcgt2=shift@line;
			my $nu="@line";
			$nu=~s/ac/ab/g;
			$nu=~s/bc/ac/g;
			$nu=~s/ /\t/g;
			print OU "$aabcid\tbc\taa\t$nu\n";
		}
	}
	unlink "parent2_abxcc.txt";
}
chdir "$INDELGT/pls";
`perl segregation_ratio.pl -v parent1.vcf -p $population -d $outDirName -i $progeny_number`;
chdir "$INDELGT/$outDirName";
if(-e "parent1_abxaa.txt"){
	`mv parent1_abxaa.txt abxaa.txt`;
}
if(-e "parent1_abxcc.txt"){
	`mv parent1_abxcc.txt abxcc.txt`;
}
my %cls;
open PC1, "parent1.cls";
while(<PC1>){
	chomp;
	$cls{$_}=1;
}
close PC1;
open PC2, "parent2.cls";
open PCOUT, ">parent2.cls1";
while(<PC2>){
	chomp;
	if(exists $cls{$_}){
		next;
	}else{
		print PCOUT "$_\n";
	}
}
close PC2;
close PCOUT;
undef %cls;
`cat parent1.cls parent2.cls1 > parent.old_cls`;
unlink "parent1.cls","parent2.cls","parent2.cls1";
open CLS, "parent.old_cls";
open CLSOUT, ">parent.cls";
while(<CLS>){
	chomp;
	my @line=split/\s+/,$_;
	my @sort=sort {(($a=~/(\d+)/)[0] <=> ($b=~/(\d+)/)[0]) or ((split/\s+/,$a)[1] <=> (split/\s+/,$b)[1])}@line;
	print CLSOUT join ("\t",@sort),"\n";
}
close CLS;
close CLSOUT;
unlink "parent.old_cls";
chdir '..';
sub prefixfq0{
        my ($fq1,$fq2) = @_;
        my $i;
        my ($fastq1,$fastq2,$fqid1,$fqid2);
        my ($n1,$n2);
        my ($s1,$s2);
        my ($prefix1,$prefix2);
        my (@tmp1,@tmp2);
        $n1 = length($fq1);
        $n2 = length($fq2);

        if($n1 != $n2){
                return (0,"");
        }

        for($i=0;$i<$n1;$i++){
                $s1 = substr($fq1,$i,1);
                $s2 = substr($fq2,$i,1);
                if($s1 ne $s2 && $s1 eq "1" && $s2 eq "2" && $i>0){
                        $prefix1 = substr($fq1,0,$i);
                        $prefix2 = substr($fq2,0,$i);
                        $fqid1 = $prefix1;
                        $fqid2 = $prefix2;
                        $fastq1=$fq1;
                        $fastq2=$fq2;
                        $prefix1 =~ s/[-\_\._]$//;
                        $prefix2 =~ s/[-\_\._]$//;
                        $fastq1 =~ s/$fqid1\s*1//;
                        $fastq2 =~ s/$fqid2\s*2//;
                        if($prefix1=~/\// or $prefix2=~/\//){
                                @tmp1 = split /\//,$prefix1;
                                @tmp2 = split /\//,$prefix2;
                                $prefix1 = pop @tmp1;
                                $prefix2 = pop @tmp2;
                        }if($prefix1 eq $prefix2 && $fastq1 eq $fastq2){
                                return (1,$prefix1);
                        }elsif($prefix1 ne $prefix2 || $fastq1 ne $fastq2){
                                return (0,"");
                        }
                }
        }

        return (0,"");
}

####################################The end of indels genotyping of the progeny#############################################
