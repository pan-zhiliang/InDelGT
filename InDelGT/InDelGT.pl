#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Parallel::ForkManager;

sub prtHelp{
	print "\nInDelGT is a software for genotyping a large number of InDel loci in hybrid populations with RAD-Seq data.\n";
	print "Contact: Chunfa Tong <tongchf\@njfu.edu.cn>\n";
	print "Version: 1.0\n";
	print "Usage: perl InDelGT.pl -o directory\n";
	print "\n";
	print "Options:\n";
	print "	-o <directory>	create a directory for storing output file of the InDel genotyping results\n";
	print "	--help|h	help\n";
}

my %prs;
my ($str1,$str2);
my $outDirName;
my $help;

if(@ARGV<1){
	prtHelp();
	exit;
}

GetOptions(
	"o:s"=>\$outDirName,
	"help|h"=>\$help
);


if($help){
	prtHelp();
	exit;
}

unless($outDirName){
	print STDERR "\nError: please create a directory for storing output file!\n\n";
	prtHelp();
	exit;
}

my $np=0;
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

my ($bwa,$samtools,$bcftools,$progeny_fold,$INDELGT,$reference,$pls,$male1,$male2,$female1,$female2,$progeny_number,$MAPQ,$threads);

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

if(exists($prs{"PROGENY_FOLD"})){
	$progeny_fold = $prs{"PROGENY_FOLD"};
}else{
	die "Error! Please check the line of PROGENY_FOLD in \"parameters.ini\".\n";
}

if(exists($prs{"INDELGT_FOLD"})){
	$INDELGT = $prs{"INDELGT_FOLD"};
}else{
	die "Error! Please check the line of INDELGT_FOLD in \"parameters.ini\".\n";
}

if(exists($prs{"REFERENCE_GENOME_FILE"})){
	$reference = $prs{"REFERENCE_GENOME_FILE"};
}else{
	die "Error! Please check the line of REFERENCE_GENOME_FILE in \"parameters.ini\".\n";
}

if(exists($prs{"PLSFOLD"})){
	$pls = $prs{"PLSFOLD"};
}else{
	die "Error! Please check the line of PLSFOLD in \"parameters.ini\".\n";
}

if(exists($prs{"MALEPARENT_1_FASTQ_FILE"})){
	$male1 = $prs{"MALEPARENT_1_FASTQ_FILE"};
}else{
	die "Error! Please check the line of MALEPARENT_1_FASTQ_FILE in \"parameters.ini\".\n";
}

if(exists($prs{"MALEPARENT_2_FASTQ_FILE"})){
	$male2 = $prs{"MALEPARENT_2_FASTQ_FILE"};
}else{
	die "Error! Please check the line of MALEPARENT_2_FASTQ_FILE in \"parameters.ini\".\n";
}

if(exists($prs{"FEMALEPARENT_1_FASTQ_FILE"})){
	$female1 = $prs{"FEMALEPARENT_1_FASTQ_FILE"};
}else{
	die "Error! Please check the line of FEMALEPARENT_1_FASTQ_FILE in \"parameters.ini\".\n";
}

if(exists($prs{"FEMALEPARENT_2_FASTQ_FILE"})){
	$female2 = $prs{"FEMALEPARENT_2_FASTQ_FILE"};
}else{
	die "Error! Please check the line of FEMALEPARENT_2_FASTQ_FILE in \"parameters.ini\".\n";
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

my @refGenome=split/\//,$reference;
`cp $reference $pls`;
chdir $pls or die "Error: can't cd to directory '$pls' : $!";
mkdir $outDirName, or die "Error: can't create directory '$outDirName' : $!" unless( -d $outDirName);
if(!-e $reference){
	print STDERR "Error: The reference genome was not prepared. Please prepare the reference genome!\n\n";
	exit;
}if(-e $refGenome[-1]){
	`$bwa index $refGenome[-1]`;
}
if((!-e $male1) or (!-e $male2)){
	print STDERR "Error: The pair-end reads of the male parent were not prepared. Please prepare the pair-end reads of the male parent!\n\n";
	exit;
}if((-e $male1) or (-e $male2)){
	`$bwa mem -t $threads $refGenome[-1] $male1 $male2 > male.sam`;
}
`$samtools sort -@ $threads -n -o male.sorted.sam male.sam`;
if(-e "male.sorted.sam"){
	unlink ("male.sam");
}
`$samtools fixmate -@ $threads -m male.sorted.sam male.sorted.fixmate`;
if(-e "male.sorted.fixmate"){
	unlink ("male.sorted.sam");
}
`$samtools sort  -@ $threads male.sorted.fixmate -o male.sorted.fixmate1`;
if(-e "male.sorted.fixmate1"){
	unlink ("male.sorted.fixmate");
}
`$samtools markdup -@ $threads  -r male.sorted.fixmate1 male.sorted.markdup`;
if(-e "male.sorted.markdup"){
	unlink ("male.sorted.fixmate1");
}
`$samtools view -@ $threads -b -S male.sorted.markdup -q $MAPQ > male.bam`;
if(-e "male.bam"){
	unlink ("male.sorted.markdup");
}
`$samtools sort -@ $threads male.bam > male.sorted.bam`;
if(-e "male.sorted.bam"){
	unlink ("male.bam");
}
`$samtools index male.sorted.bam`;
`$bcftools mpileup -Obuzv -a  AD,INFO/AD -f $refGenome[-1] male.sorted.bam > male.bcf`;
`$bcftools call -m -v -f gq male.bcf > male.vcf`;

####################################Identification of female parent indels###############################################

if((!-e $female1) or (!-e $female2)){
	print STDERR "Error: The pair-end reads of the female parent were not prepared. Please prepare the pair-end reads of the female parent!\n\n";
	exit;
}if((-e $female1) or (-e $female2)){
	`$bwa mem -t $threads $refGenome[-1] $female1 $female2 > female.sam`;
}
`$samtools sort -@ $threads -n -o female.sorted.sam female.sam`;
if(-e "female.sorted.sam"){
	unlink ("female.sam");
}
`$samtools fixmate -@ $threads -m female.sorted.sam female.sorted.fixmate`;
if(-e "female.sorted.fixmate"){
	unlink ("female.sorted.sam");
}
`$samtools sort  -@ $threads female.sorted.fixmate -o female.sorted.fixmate1`;
if(-e "female.sorted.fixmate1"){
	unlink ("female.sorted.fixmate");
}
`$samtools markdup -@ $threads  -r female.sorted.fixmate1 female.sorted.markdup`;
if(-e "female.sorted.markdup"){
	unlink ("female.sorted.fixmate1");
}
`$samtools view -@ $threads -b -S female.sorted.markdup -q $MAPQ > female.bam`;
if(-e "female.bam"){
	unlink ("female.sorted.markdup");
}
`$samtools sort -@ $threads female.bam > female.sorted.bam`;
if(-e "female.sorted.bam"){
        unlink ("female.bam");
}
`$samtools index female.sorted.bam`;
`$bcftools mpileup -Obuzv -a  AD,INFO/AD -f $refGenome[-1] female.sorted.bam > female.bcf`;
`$bcftools call -m -v -f gq female.bcf > female.vcf`;

####################################The end of identification of parents indels#############################################

##########################################Indels genotyping of the parents##################################################

`perl $pls\/parent_genotyping.pl -v male.vcf -b female.bcf -o $outDirName`;
`perl $pls\/parent_genotyping.pl -v female.vcf -b male.bcf -o $outDirName`;
`perl $pls\/homozygote_filter.pl -v male_female_loc.parentgt -b male.sorted.bam -o $outDirName`;
`perl $pls\/homozygote_filter.pl -v female_male_loc.parentgt -b female.sorted.bam -o $outDirName`;

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
		$str1 = "$progeny_fold\/$prg1fq[$i-1]";
		$str2 = "$progeny_fold\/$prg2fq[$i-1]";
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
				`$bwa mem -t $threads $refGenome[-1] $progeny_fold/$fq[0] $progeny_fold/$fq[1] > $fqprfx[$i-1]_progeny.sam`;
				`$samtools sort -@ $threads -n -o $fqprfx[$i-1]_progeny.sorted.sam $fqprfx[$i-1]_progeny.sam`;
				if(-e "$fqprfx[$i-1]_progeny.sorted.sam"){
					unlink ("$fqprfx[$i-1]_progeny.sam");
				}
				`$samtools fixmate -@ $threads -m $fqprfx[$i-1]_progeny.sorted.sam $fqprfx[$i-1]_progeny.sorted.fixmate`;
				if(-e "$fqprfx[$i-1]_progeny.sorted.fixmate"){
					unlink ("$fqprfx[$i-1]_progeny.sorted.sam");
				}
				`$samtools sort  -@ $threads -o $fqprfx[$i-1]_progeny.sorted.fixmate1 $fqprfx[$i-1]_progeny.sorted.fixmate`;
				if(-e "$fqprfx[$i-1]_progeny.sorted.fixmate1"){
					unlink ("$fqprfx[$i-1]_progeny.sorted.fixmate");
				}
				`$samtools markdup -@ $threads  -r $fqprfx[$i-1]_progeny.sorted.fixmate1 $fqprfx[$i-1]_progeny.sorted.markdup`;
				if(-e "$fqprfx[$i-1]_progeny.sorted.markdup"){
					unlink ("$fqprfx[$i-1]_progeny.sorted.fixmate1");
				}
				`$samtools view -@ $threads -b -S $fqprfx[$i-1]_progeny.sorted.markdup -q $MAPQ > $fqprfx[$i-1]_progeny.bam`;
				if(-e "$fqprfx[$i-1]_progeny.bam"){
				unlink ("$fqprfx[$i-1]_progeny.sorted.markdup");
				}
				`$samtools sort -@ $threads $fqprfx[$i-1]_progeny.bam > $fqprfx[$i-1]_progeny.sorted.bam`;
				if(-e "$fqprfx[$i-1]_progeny.sorted.bam"){
					unlink ("$fqprfx[$i-1]_progeny.bam");
				}
				`$samtools index $fqprfx[$i-1]_progeny.sorted.bam`;
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
	`$bcftools mpileup -Obuzv -a  AD,INFO/AD -f $refGenome[-1] $Sorted_Name\.sorted.bam > $Sorted_Name.bcf`;
	`perl $pls\/progeny_genotyping.pl -v male.vcf -b $Sorted_Name\.bcf -o $outDirName`;
	`perl $pls\/progeny_genotyping.pl -v female.vcf -b $Sorted_Name\.bcf -o $outDirName`;
	`perl $pls\/homozygote_filter.pl -v $Sorted_Name\_male_loc.gt -b $Sorted_Name\.sorted.bam -o $outDirName`;
	`perl $pls\/homozygote_filter.pl -v $Sorted_Name\_female_loc.gt -b $Sorted_Name\.sorted.bam -o $outDirName`;
	$pm->finish;
}
$pm->wait_all_children;

unlink "$refGenome[-1]";
unlink glob"$refGenome[-1].*";
if(-e "nohup.out"){
	`mv nohup.out $outDirName`;
}
`mv *.bcf $outDirName`;
`mv *.sorted.bam $outDirName`;
`mv *.sorted.bam.bai $outDirName`;
`mv *.vcf $outDirName`;
`mv *.bcf.csi $outDirName`;
`mv *.cls $outDirName`;
`perl $pls\/segregation_ratio.pl -v male.vcf -d $outDirName -i $progeny_number`;
chdir $outDirName;
`mv male_abxaa.txt aaxab.txt`;
if(-e "male_abxcc.txt"){
	open IN, "male_abxcc.txt";
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
	unlink "male_abxcc.txt";
}
chdir '..';
`perl $pls\/segregation_ratio.pl -v female.vcf -d $outDirName -i $progeny_number`;
chdir $outDirName;
`mv female_abxaa.txt abxaa.txt`;
if(-e "female_abxcc.txt"){
	`mv female_abxcc.txt abxcc.txt`;
}
chdir '..';
chdir '..';
if(-d $outDirName){
	`rm -r $outDirName`;
	`mv $pls/$outDirName ./`;
}if(!-d $outDirName){
	`mv $pls/$outDirName ./`;
}

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
