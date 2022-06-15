#!/usr/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;
use Parallel::ForkManager;

sub prtHelp{
	print "\nInDelGT is an integrated pipeline for extracting InDel genotypes in a hybrid population with next generation DNA sequencing data.\n";
	print "Contact: Chunfa Tong <tongchf\@njfu.edu.cn>\n";
	print "Version: 1.0\n";
	print "Usage: perl InDelGT.pl [option] <type> -o <directory>\n";
	print "\n";
	print "Options:\n";
	print "	--nocatalog	skip the step for generating parental InDel catalogs\n";
	print "	--nocall	skip the step for calling InDel genotypes for all progeny\n";
	print "	--nofilter	skip the step for analyzing InDel segregation types\n";
	print "	-p  <type>	the type of population: CP;BC1;BC2;F2 (default: CP)\n";
	print "	-o  <dir>	create a directory for storing output file of the InDel genotyping results\n";
	print "	--help|h	help\n";
}

my %prs;
my %hash;
my ($str1,$str2);
my $outDirName;
my $population='CP';
my $catalog=1;
my $call=1;
my $filter=1;
my $help;
my $np=0;
if(@ARGV<1){
	prtHelp();
	exit;
}

GetOptions(
	"p:s"=>\$population,
	"o:s"=>\$outDirName,
	"catalog!"=>\$catalog,
	"call!"=>\$call,
	"filter!"=>\$filter,
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

my $cwd=getcwd;
open(PRS,"$cwd/parameters.ini") || die "\nError\! This program needs a parameter file,namely \"parameters.ini\" in $cwd.\n\n";
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

my $i;
my $str0;
my $pm;
my %cls;
my ($flag,$flag1);
my @prg1fq;
my @prg2fq;
my $parent_threads;
my ($fqprfx0,@fqprfx);
my ($fqprfx1,@fqprfxparent);
my ($bwa,$samtools,$bcftools,$INDELGT,$datafold,$reference,$parent1,$parent2,$progeny_number,$MAPQ,$threads,$HOMOZYGOUS,$HETEROZYGOUS,$GQ,$MISPCT,$PVALUE);

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
        unless( $MAPQ =~ /\d+/ ){
                die "Error\! The number of MAPQ is wrong.\n";
        }
}else{
        die "Error! Please check the line of MAPQ in \"parameters.ini\".\n";
}

if(exists($prs{"THREADS"})){
        $threads = $prs{"THREADS"};
        unless( $threads =~ /\d+/ ){
                die "Error\! The number of threads is wrong.\n";
        }
}else{
        die "Error! Please check the line of THREADS in \"parameters.ini\".\n";
}

if(exists($prs{"HOMOZYGOUS_DEPTH"})){
        $HOMOZYGOUS = $prs{"HOMOZYGOUS_DEPTH"};
        unless( $HOMOZYGOUS =~ /\d+/ ){
                die "Error\! The value of homozygous depth is wrong.\n";
        }
}else{
        die "Error! Please check the line of HOMOZYGOUS_DEPTH in \"parameters.ini\".\n";
}

if(exists($prs{"HETEROZYGOUS_DEPTH"})){
        $HETEROZYGOUS = $prs{"HETEROZYGOUS_DEPTH"};
        unless( $HETEROZYGOUS =~ /\d+/ ){
                die "Error\! The value of heterozygous depth is wrong.\n";
        }
}else{
        die "Error! Please check the line of HETEROZYGOUS_DEPTH in \"parameters.ini\".\n";
}

if(exists($prs{"GQ"})){
        $GQ = $prs{"GQ"};
        unless( $GQ =~ /\d+/ ){
                die "Error\! The value of genotype quality (GQ) is wrong.\n";
        }
}else{
        die "Error! Please check the line of GQ in \"parameters.ini\".\n";
}

if(exists($prs{"MISPCT"})){
        $MISPCT = $prs{"MISPCT"};
        unless($MISPCT =~ /\d+/){
                die "Error\! The number for the percent of missing genotypes is wrong.\n";
        }
        if($MISPCT < 0 || $MISPCT > 100){
                die "Error\! The parameter of MISPCT is wrong.\n";
        }
}else{
        die "Error! Please check the line of MISPCT in \"parameters.ini\".\n";
}

if(exists($prs{"PVALUE"})){
        $PVALUE = $prs{"PVALUE"};
        if($PVALUE < 0 || $PVALUE > 1.0){
                die "Error\! The parameter of PVALUE  is wrong.\n";
        }
}else{
        die "Error! Please check the line of PVALUE in \"parameters.ini\".\n";
}

####################################Identification of male parent indels###############################################
my $home=(`echo \$HOME`);
chomp$home;
$bwa=~s/^\~/$home/g if $bwa=~/^\~/;
$samtools=~s/^\~/$home/g if $samtools=~/^\~/;
$bcftools=~s/^\~/$home/g if $bcftools=~/^\~/;
$datafold=~s/^\~/$home/g if $datafold=~/^\~/;
$INDELGT=~s/^\~/$home/g if $INDELGT=~/^\~/;
$bwa=~s/\/$//g if $bwa=~/\/$/;
$samtools=~s/\/$//g if $samtools=~/\/$/;
$bcftools=~s/\/$//g if $bcftools=~/\/$/;
$INDELGT=~s/\/$//g if $INDELGT=~/\/$/;
$datafold=~s/\/$//g if $datafold=~/\/$/;
die "Error! Please check the line of BWA_FOLD in \"parameters.ini\".\n" if (!-d $bwa);
die "Error! Please check the line of SAMTOOLS_FOLD in \"parameters.ini\".\n" if (!-d $samtools);
die "Error! Please check the line of BCFTOOLS_FOLD in \"parameters.ini\".\n" if (!-d $bcftools);
die "Error! Please check the line of DATAFOLD in \"parameters.ini\".\n" if (!-d $datafold);
die "Error! Please check the line of INDELGT_FOLD in \"parameters.ini\".\n" if (!-d $INDELGT);
my $ref="$datafold/$reference";
my @parent1_fq=split/\s+/,$parent1;
my @parent2_fq=split/\s+/,$parent2;
my @refGenome=split/\//,$ref;
my $parent1_fq_file1="$datafold/$parent1_fq[0]";
my $parent1_fq_file2="$datafold/$parent1_fq[1]";
my $parent2_fq_file1="$datafold/$parent2_fq[0]";
my $parent2_fq_file2="$datafold/$parent2_fq[1]";
if($outDirName=~/\/$/){
	$outDirName=~s/\/$//g;
}
$outDirName=~s/^\~/$home/g if $outDirName =~/^\~/;
if($outDirName=~/^\./ or $outDirName=~/^\.\./ or $outDirName!~/\//){
        $outDirName="$cwd/$outDirName";
}
mkdir $outDirName, or die "Error: can't create directory '$outDirName' : $!" unless( -d $outDirName);
chdir $outDirName or die "Error: can't cd to directory '$outDirName' : $!";
if(!-e $ref){
	print STDERR "Error: The reference genome was not prepared. Please prepare the reference genome!\n\n";
	exit;
}
if((!-e $parent1_fq_file1) or (!-e $parent1_fq_file2)){
	print STDERR "Error: The pair-end reads of the parent1 were not prepared. Please prepare the pair-end reads of the parent1!\n\n";
	exit;
}
if((!-e $parent2_fq_file1) or (!-e $parent2_fq_file2)){
	print STDERR "Error: The pair-end reads of the parent2 were not prepared. Please prepare the pair-end reads of the parent2!\n\n";
	exit;
}
($flag1,$fqprfx1) = &prefixfq0($parent1_fq[0],$parent1_fq[1]);
if($flag1 == 0){
	die "Error: the two fastq file names of the parent1 are not consistent!\n";
}
if($flag1 == 1){
	push @fqprfxparent,$fqprfx1;
}
($flag1,$fqprfx1) = &prefixfq0($parent2_fq[0],$parent2_fq[1]);
if($flag1 == 0){
        die "Error: the two fastq file names of the parent2 are not consistent!\n";
}
if($flag1 == 1){
        push @fqprfxparent,$fqprfx1;
}
my @parent_fq=([$parent1_fq_file1,$parent1_fq_file2],[$parent2_fq_file1,$parent2_fq_file2]);
if($catalog==1){
	if(!-e "$INDELGT/pls/parent_genotyping.pl"){
		die "Error: not find the Perl program parent_genotyping.pl in $INDELGT/pls.\n";
	}
	`cp $ref ./`;
	if($bwa=~/bwa-mem2/){
		`$bwa/bwa-mem2 index $reference`;
	}else{
		`$bwa/bwa index $reference`;
	}
	$parent_threads="$threads";
	if($threads>2){
		$parent_threads=2;
	}
	$pm = new Parallel::ForkManager($parent_threads);	
	for($i=0;$i<2;$i++){
		my $pid = $pm->start and next;
		if($bwa=~/bwa-mem2/){
			`$bwa/bwa-mem2 mem $reference $parent_fq[$i][0] $parent_fq[$i][1] > $fqprfxparent[$i].sam`;
		}else{
			`$bwa/bwa mem $reference $parent_fq[$i][0] $parent_fq[$i][1] > $fqprfxparent[$i].sam`;
		}
		`$samtools/samtools sort -n -o $fqprfxparent[$i].sorted.sam $fqprfxparent[$i].sam`;
		if(-e "$fqprfxparent[$i].sorted.sam"){
			unlink ("$fqprfxparent[$i].sam");
		}
		`$samtools/samtools fixmate -m $fqprfxparent[$i].sorted.sam $fqprfxparent[$i].sorted.fixmate`;
		if(-e "$fqprfxparent[$i].sorted.fixmate"){
			unlink ("$fqprfxparent[$i].sorted.sam");
		}
		`$samtools/samtools sort $fqprfxparent[$i].sorted.fixmate -o $fqprfxparent[$i].sorted.fixmate1`;
		if(-e "$fqprfxparent[$i].sorted.fixmate1"){
			unlink ("$fqprfxparent[$i].sorted.fixmate");
		}
		`$samtools/samtools markdup -r $fqprfxparent[$i].sorted.fixmate1 $fqprfxparent[$i].sorted.markdup`;
		if(-e "$fqprfxparent[$i].sorted.markdup"){
			unlink ("$fqprfxparent[$i].sorted.fixmate1");
		}
		`$samtools/samtools view -b -S $fqprfxparent[$i].sorted.markdup -q $MAPQ > $fqprfxparent[$i].bam`;
		if(-e "$fqprfxparent[$i].bam"){
			unlink ("$fqprfxparent[$i].sorted.markdup");
		}
		`$samtools/samtools sort $fqprfxparent[$i].bam > $fqprfxparent[$i].sorted.bam`;
		if(-e "$fqprfxparent[$i].sorted.bam"){
			unlink ("$fqprfxparent[$i].bam");
		}
		`$samtools/samtools index $fqprfxparent[$i].sorted.bam`;
		`$bcftools/bcftools mpileup -Obuzv -a  AD,INFO/AD -f $reference $fqprfxparent[$i].sorted.bam > $fqprfxparent[$i].bcf`;
		`$bcftools/bcftools call -m -v -f gq $fqprfxparent[$i].bcf > $fqprfxparent[$i].vcf`;
		$pm->finish;
	}
	$pm->wait_all_children;

####################################The end of identification of parents indels#############################################

##########################################Indels genotyping of the parents##################################################

	`perl $INDELGT/pls/parent_genotyping.pl -q $GQ -c $bcftools -ho $HOMOZYGOUS -he $HETEROZYGOUS -v $fqprfxparent[0].vcf -b $fqprfxparent[1].bcf -o $outDirName`;
	`perl $INDELGT/pls/parent_genotyping.pl -q $GQ -c $bcftools -ho $HOMOZYGOUS -he $HETEROZYGOUS -v $fqprfxparent[1].vcf -b $fqprfxparent[0].bcf -o $outDirName`;
	for ($i=0;$i<2;$i++){
		open IN, "$fqprfxparent[$i].vcf";
		open OU, ">01-$fqprfxparent[$i].vcf";
		while(<IN>){
        		chomp;
		        if($_=~/#/ or $_!~/INDEL/){
        		        next;
	        	}elsif($_!~/#/){
        	        	my @line=split/\s+/,$_;
	                	my @vcf_gt=split/\:/,$line[-1];
		                my @vcf_allele=split/\//,$vcf_gt[0];
        		        if($vcf_allele[0]!=$vcf_allele[1]){
                		        print OU "$_\n";
		                }else{
        		                next;
	        	        }
	        	}else{
        	        	next;
        		}
		}
		close IN;
		close OU;
		open IN,"01-$fqprfxparent[$i].vcf";
		open OU,">01-$fqprfxparent[$i].vcf.gt";
		while(<IN>){
		        chomp;
		        my @line=split/\s+/,$_;
			my $len1=(length$line[3]);
			my $len2=($line[1]+$len1-1);
		        my @b=split/\:/,$line[-1];
		        my @allel=split/,/,$line[4];
		        my @dp=split /;/,$line[7];
		        my @cd=split /:/,$line[-1];
		        my @dp1=split /=/,$dp[3];
		        my @dp2=split /=/,$dp[4];
		        my @dp3=split /,/,$dp2[1];
		        $cls{0}=$line[3];
		        my $i=1;
		        foreach (@allel){
		                $cls{$i}=$_;
		                $i++;
		        }
		        my$GT=$b[0];
		        if(($GT!~/0\/0/ && $GT!~/\.\/\./) && (($dp1[1]<$HOMOZYGOUS) or ($b[-1]<$GQ) or ($line[5]<20))){
		                print OU "$line[0]:$line[1]-$len2\t$line[3]\t$line[4]\t_\/_\n";
		        }elsif($GT=~/0\/0/ && $line[5]>=20){
		                print OU "$line[0]:$line[1]-$len2\t$line[3]\t$line[4]\t$GT\n";
		        }elsif($GT=~/\.\/\./){
		                print OU "$line[0]:$line[1]-$len2\t$line[3]\t$line[4]\t_/_\n";
		        }elsif($GT=~/1\/1/ && (($dp3[1]<$HOMOZYGOUS) or ($b[-1]<$GQ) or $line[5]<20)){
		                print OU "$line[0]:$line[1]-$len2\t$line[3]\t$line[4]\t_/_\n";
		        }elsif($GT=~/0\/1/ && (($dp3[0]<$HETEROZYGOUS) or ($dp3[1]<$HETEROZYGOUS) or ($b[-1]<$GQ) or $line[5]<20)){
		                print OU "$line[0]:$line[1]-$len2\t$line[3]\t$line[4]\t_/_\n";
		        }elsif($GT!~/0\/0/){
		                $b[0]=~s:(\d).*(\d):$cls{$1}/$cls{$2}:g;
		                print OU "$line[0]:$line[1]-$len2\t$line[3]\t$line[4]\t$b[0]\n";
        		}
		        undef %cls;
		}
		close IN;
		close OU;
	}
	`cat 01-$fqprfxparent[0].vcf.gt $fqprfxparent[0].gt |sort |uniq > $fqprfxparent[0].gt1`;
	`cat 01-$fqprfxparent[1].vcf.gt $fqprfxparent[1].gt |sort |uniq > $fqprfxparent[1].gt1`;
	unlink "01-$fqprfxparent[0].vcf","01-$fqprfxparent[1].vcf","01-$fqprfxparent[0].vcf.gt","01-$fqprfxparent[1].vcf.gt","$fqprfxparent[0].gt","$fqprfxparent[1].gt";
	`mv $fqprfxparent[0].gt1 $fqprfxparent[0].gt`;
	`mv $fqprfxparent[1].gt1 $fqprfxparent[1].gt`;
	`cat $fqprfxparent[0]-$fqprfxparent[1].vcf.loci-seq $fqprfxparent[1]-$fqprfxparent[0].vcf.loci-seq |sort |uniq > parent_total.seq`;
	open IN,"parent_total.seq";
	open OU,">parent_total.seq.indeluniq";
	while(<IN>){
		my @line=split/\s+/,$_;
		my $pos_id="$line[0]\t$line[1]";
		next if exists $hash{$pos_id};
		$hash{$pos_id}= 1;
		print OU "$_";
	}
	close IN;
	close OU;
	undef %hash;
	`sort parent_total.seq parent_total.seq.indeluniq |uniq -u > uniq.loc`;
	open IN,"uniq.loc";
	while(<IN>){
		chomp;
		my @line=split/\s+/,$_;
		$hash{"$line[0]\t$line[1]"}=1;
	}
	open IN, "parent_total.seq";
	open OU, ">parent.cls";
	while(<IN>){
		chomp;
		my @line=split/\s+/,$_;
		if(exists $hash{"$line[0]\t$line[1]"} ){
			next;
		}else{
			print OU "$_\n";
		}
	}
	close IN;
	close OU;
	undef %hash;
	unlink "$fqprfxparent[0]-$fqprfxparent[1].vcf.loci-seq","$fqprfxparent[1]-$fqprfxparent[0].vcf.loci-seq","uniq.loc","parent_total.seq","parent_total.seq.indeluniq";
			
}
####################################The end of indels genotyping of the parents#############################################

##########################################Indels genotyping of the progeny##################################################
my $profq;
my @prbamname;
if($call==1){
	if(!-e glob ("$reference*")){
		die "Error: The reference index file was not generated. Please perform the 'catalog' step!\n\n";
	}
	if(!-e "parent.cls"){
		die "Error: The parent catalog file parent.cls was not generated. Please perform the 'catalog' step!\n\n";
	}
	if(!-e "$INDELGT/pls/progeny_genotyping.pl"){
		die "Error: not find the Perl program progeny_genotyping.pl in $INDELGT/pls.\n";
	}
	if(!-e "$fqprfxparent[0].vcf"){
		die "Error: The parent1 VCF file $fqprfxparent[0].vcf was not generated. Please perform the 'catalog' step!\n\n";
	}
	if(!-e "$fqprfxparent[1].vcf"){
                die "Error: The parent2 VCF file $fqprfxparent[1].vcf was not generated. Please perform the 'catalog' step!\n\n";
        }
	if(!-e "$fqprfxparent[0].gt"){
		die "Error: The parent1 InDel genotypes file $fqprfxparent[0].gt was not generated. Please perform the 'catalog' step!\n\n";
	}
        if(!-e "$fqprfxparent[1].gt"){
                die "Error: The parent2 InDel genotypes file $fqprfxparent[1].gt was not generated. Please perform the 'catalog' step!\n\n";
        }
	my $progeny_threads=$threads;
	if($progeny_number<$threads){
		$progeny_threads=$progeny_number;
	}
	$pm = new Parallel::ForkManager($progeny_threads);
	for($i = 1;$i <= $np;$i++){
		my $pid = $pm->start and next;
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
			}
		}else{
			die "Error! Please check the line of PROGENY$i in \"parameters.ini\".\n";
		}
		my @fq=split/\s+/,$prs{$str0};
		($flag,$profq) = &prefixfq0($fq[0],$fq[1]);
		if($bwa=~/bwa-mem2/){
			`$bwa/bwa-mem2 mem $reference $datafold/$fq[0] $datafold/$fq[1]  > $profq.sam`;
		}else{
			`$bwa/bwa mem $reference $datafold/$fq[0] $datafold/$fq[1] > $profq.sam`;
		}
		`$samtools/samtools sort -n -o $profq.sorted.sam $profq.sam`;
		if(-e "$profq.sorted.sam"){
			unlink ("$profq.sam");
		}
		`$samtools/samtools fixmate -m $profq.sorted.sam $profq.sorted.fixmate`;
		if(-e "$profq.sorted.fixmate"){
			unlink ("$profq.sorted.sam");
		}
		`$samtools/samtools sort -o $profq.sorted.fixmate1 $profq.sorted.fixmate`;
		if(-e "$profq.sorted.fixmate1"){
			unlink ("$profq.sorted.fixmate");
		}
		`$samtools/samtools markdup -r $profq.sorted.fixmate1 $profq.sorted.markdup`;
		if(-e "$profq.sorted.markdup"){
			unlink ("$profq.sorted.fixmate1");
		}
		`$samtools/samtools view -b -S $profq.sorted.markdup -q $MAPQ > $profq.bam`;
		if(-e "$profq.bam"){
			unlink ("$profq.sorted.markdup");
		}
		`$samtools/samtools sort $profq.bam > $profq.sorted.bam`;
		if(-e "$profq.sorted.bam"){
			unlink ("$profq.bam");
		}
		`$samtools/samtools index $profq.sorted.bam`;
		`$bcftools/bcftools mpileup -Obuzv -a  AD,INFO/AD -f $reference $profq.sorted.bam > $profq.bcf`;
		`perl $INDELGT/pls/progeny_genotyping.pl -q $GQ -c $bcftools -ho $HOMOZYGOUS -he $HETEROZYGOUS -v parent.cls -b $profq\.bcf -o $outDirName`;
		$pm->finish;
	}
	$pm->wait_all_children;
}
if($filter==1){
        if(!-e "$INDELGT/pls/segregation_ratio.pl"){
                die "Error: not find the Perl program segregation_ratio.pl in $INDELGT/pls.\n";
        }
	if(!-e "$fqprfxparent[0].vcf"){
		die "Error: The parent1 VCF file $fqprfxparent[0].vcf was not generated. Please perform the 'catalog' step!\n\n";
	}
	if(!-e "$fqprfxparent[1].vcf"){
		die "Error: The parent2 VCF file $fqprfxparent[1].vcf was not generated. Please perform the 'catalog' step!\n\n";
	}
        if(!-e "$fqprfxparent[0].gt"){
                die "Error: The parent1 InDel genotypes file $fqprfxparent[0].gt was not generated. Please perform the 'catalog' step!\n\n";
        }
        if(!-e "$fqprfxparent[1].gt"){
                die "Error: The parent2 InDel genotypes file $fqprfxparent[1].gt was not generated. Please perform the 'catalog' step!\n\n";
        }
        for($i = 1;$i <= $np;$i++){
                $str0 = "PROGENY$i";
                if(exists($prs{$str0})){
                        ($prg1fq[$i-1],$prg2fq[$i-1]) = split /\s+/,$prs{$str0};
                        ($flag,$fqprfx0) = &prefixfq0($prg1fq[$i-1],$prg2fq[$i-1]);
                        if($flag == 1){
                                push @fqprfx,$fqprfx0;
                        }else{
                                die "Error: The file names of $prg1fq[$i-1] and $prg2fq[$i-1] are not consistent!\n";
                                exit;
                        }
		}
		if(!-e "$fqprfx[$i-1].gt"){
			die "Error: The progeny InDel genotypes file $fqprfx[$i-1].gt was not generated. Please perform the calling step!\n\n";
		}
	}
	if(!-e "parameters.ini"){
		`mv $cwd/parameters.ini ./`;
	}
        if($population eq 'BC1'){
                unlink "aaxab.txt" if (-e "aaxab.txt");
                unlink "abxab.txt" if (-e "abxab.txt");
                unlink "abxcc.txt" if (-e "abxcc.txt");
                unlink "aaxbc.txt" if (-e "aaxbc.txt");
                unlink "abxac.txt" if (-e "abxac.txt");
                unlink "abxcd.txt" if (-e "abxcd.txt");
        }
        if($population eq 'BC2'){
                unlink "abxaa.txt" if (-e "abxaa.txt");
                unlink "abxab.txt" if (-e "abxab.txt");
                unlink "abxcc.txt" if (-e "abxcc.txt");
                unlink "aaxbc.txt" if (-e "aaxbc.txt");
                unlink "abxac.txt" if (-e "abxac.txt");
                unlink "abxcd.txt" if (-e "abxcd.txt");
        }
        if($population eq 'F2'){
                unlink "abxaa.txt" if (-e "abxaa.txt");
                unlink "aaxab.txt" if (-e "aaxab.txt");
                unlink "abxcc.txt" if (-e "abxcc.txt");
                unlink "aaxbc.txt" if (-e "aaxbc.txt");
                unlink "abxac.txt" if (-e "abxac.txt");
                unlink "abxcd.txt" if (-e "abxcd.txt");
        }
	unless($population eq "BC1"){
		`perl $INDELGT/pls/segregation_ratio.pl -v $fqprfxparent[1].vcf -p $population -o $outDirName -i $progeny_number -q $GQ -ho $HOMOZYGOUS -he $HETEROZYGOUS -m $MISPCT -a $PVALUE`;
	}
	if(-e "$fqprfxparent[1]_abxaa.txt"){
		open IN, "$fqprfxparent[1]_abxaa.txt";
		open OU, ">aaxab.txt";
		while(<IN>){
			chomp;
			my @line=split/\s+/,$_;
			my $pos=shift@line;
			my $p1=shift@line;
			my $p2=shift@line;
			my $p3=shift@line;
			my $p4=shift@line;
			my $p5=shift@line;
			my $p6=shift@line;
			my $p7=shift@line;
			my $p8=shift@line;
			my $p9=shift@line;
			print OU join ("\t",$pos,$p1,$p2,$p3,$p4,$p5,$p6,$p7,$p9,$p8,@line),"\n";
		}
		close IN;
		close OU;
		unlink "$fqprfxparent[1]_abxaa.txt";
	}
	if(-e "$fqprfxparent[1]_abxcc.txt"){
		open IN, "$fqprfxparent[1]_abxcc.txt";
		open OU, ">aaxbc.txt";
		while(<IN>){
			chomp;
			if($_=~/Position/){
				my @line=split/\s+/,$_;
				my $p1=shift@line;
				my $p2=shift@line;
				my $p3=shift@line;
				my $p4=shift@line;
				my $p5=shift@line;
				my $p6=shift@line;
				my $p7=shift@line;
				my $p8=shift@line;
				my $p9=shift@line;
				my $p10=shift@line;
				my $p11=shift@line;
				print OU join ("\t",$p1,$p2,$p3,$p4,"ab","ac",$p7,$p8,$p9,$p11,$p10,@line),"\n";
			}else{
				my @line=split/\s+/,$_;
				my $aabcgt1=shift@line;
				my $aabcgt2=shift@line;
				my $aabcgt3=shift@line;
				my $aabcgt4=shift@line;
				my $aabcgt5=shift@line;
				my $aabcgt6=shift@line;
				my $aabcgt7=shift@line;
				my $aabcgt8=shift@line;
				my $aabcgt9=shift@line;
				my $aabcgt10=shift@line;
				my $aabcgt11=shift@line;
				my $nu="@line";
				$nu=~s/ac/ab/g;
				$nu=~s/bc/ac/g;
				$nu=~s/ /\t/g;
				print OU join ("\t",$aabcgt1,$aabcgt4,$aabcgt2,$aabcgt3,$aabcgt5,$aabcgt6,$aabcgt7,$aabcgt8,$aabcgt9,"aa","bc",$nu),"\n";
			}
		}
		close IN;
		close OU;
		unlink "$fqprfxparent[1]_abxcc.txt";
	}
	unless($population eq "BC2"){
		`perl $INDELGT/pls/segregation_ratio.pl -v $fqprfxparent[0].vcf -p $population -o $outDirName -i $progeny_number -q $GQ -ho $HOMOZYGOUS -he $HETEROZYGOUS -m $MISPCT -a $PVALUE`;
	}
	if(-e "$fqprfxparent[0]_abxaa.txt"){
		`mv $fqprfxparent[0]_abxaa.txt abxaa.txt`;
	}
	if(-e "$fqprfxparent[0]_abxcc.txt"){
		`mv $fqprfxparent[0]_abxcc.txt abxcc.txt`;
	}
	chdir "$cwd";
	if(!-e "parameters.ini"){
		`mv $outDirName/parameters.ini $cwd`;
	}
}
chdir "$cwd";
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
                        $prefix1 =~ s/[-\_\._a-zA-Z]$//;
                        $prefix2 =~ s/[-\_\._a-zA-Z]$//;
			if($prefix1=~/[\.\_\-\._]$/){
				$prefix1 =~ s/[\.\_\-\._]$//;
			}
			if($prefix2 =~ /[\.\_\-\._]$/){
				$prefix2 =~ s/[\.\_\-\._]$//;
			}
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
