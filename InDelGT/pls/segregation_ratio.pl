#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Statistics::Distributions;
sub prtHelp{
	print "\nThis program is used for Mendelian separation ratio analysis of indel loci\n";
	print "Contact: Chunfa Tong  <tongchf\@njfu.edu.cn>\n";
	print "Usage: perl segregation-ratio.pl [Options] <*.vcf> [Options] <directory> [Options] <int>\n";
	print "\n";
	print "Options:\n";
	print "	-v  <file>	the VCF input file for one of the two parents\n";
	print "	-i  <int>       the number of progeny\n";
	print "	-p  <type>      the type of population:CP;BC1;BC2;F2 (default: CP)\n";
	print " -q  <int>	the minimum score of InDel genotypes\n";
	print " -ho  <int>	the depth of homozygous\n";
	print " -he  <int>	the depth of heterozygous\n";
	print " -m  <int>	the percent of missing genotypes at each InDel\n";
	print " -a  <int>	the P-value\n";
	print "	-o  <directory>	the directory for storing output file of the genotyping results\n";
	print "	--help|h	help\n";

}

my $inputfile;
my $directory;
my $progeny_number;
my $population='CP';
my $help;
my %prs;
my $GQ;
my $HETEROZYGOUS;
my $HOMOZYGOUS;
my $INDELGT;
my $MISPCT;
my $PVALUE;

if(@ARGV<8){
	prtHelp();
	exit;
}

GetOptions(
	"q:s"=>\$GQ,
	"ho:s"=>\$HOMOZYGOUS,
	"he:s"=>\$HETEROZYGOUS,
	"m:s"=>\$MISPCT,
	"a:s"=>\$PVALUE,
	"v:s"=>\$inputfile,
	"i:s"=>\$progeny_number,
	"p:s"=>\$population,
	"o:s"=>\$directory,
	"help|h" => \$help
);

if ($help){
	prtHelp();
	exit;
}

unless($inputfile){
	print STDERR "\nError: the VCF file of the parent was not provided!\n\n";
	prtHelp();
	exit;
}

unless($directory){
	print STDERR "\nError: the directory for storing output file of the genotyping results was not provided!\n\n";
	prtHelp();
	exit;
}

unless($progeny_number){
	print STDERR "\nError: the number of progeny was not provided!\n\n";
	prtHelp();
	exit;
}

unless($population eq 'CP' or $population eq 'BC1' or $population eq 'BC2' or $population eq 'F2'){
	print STDERR "\nError: option '-p' only be set to 'CP' or 'BC1' or 'BC2' or 'F2'!\n\n";
	prtHelp();
	exit;
}

unless($GQ){
        print STDERR "\nError: the minimum score of InDel genotypes was not provided!\n\n";
        prtHelp();
        exit;
}

unless($HOMOZYGOUS){
        print STDERR "\nError: the depth of homozygate was not provided!\n\n";
        prtHelp();
        exit;
}

unless($HETEROZYGOUS){
        print STDERR "\nError: the depth of heterozygate was not provided!\n\n";
        prtHelp();
        exit;
}

unless($MISPCT){
        print STDERR "\nError: the percent of missing genotypes at each InDel was not provided!\n\n";
        prtHelp();
        exit;
}

unless($PVALUE){
        print STDERR "\nError: the P-value was not provided!\n\n";
        prtHelp();
        exit;
}

mkdir $directory or die "Error: can't create directory '$directory' : $!" unless(-d $directory);
chdir $directory or die "Error: can't cd to directory '$directory' : $!";
$inputfile=~s/\.vcf//g;
my %genotype;
open IN, "$inputfile.vcf";
open OU, ">01-$inputfile.vcf";
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
open IN,"01-$inputfile.vcf";
open OU,">01-$inputfile.vcf.gt";
while(<IN>){
        chomp;
        my @line=split/\s+/,$_;
        my @b=split/\:/,$line[-1];
	my @allel=split/,/,$line[4];
	my @dp=split /;/,$line[7];
	my @cd=split /:/,$line[-1];
	my @dp1=split /=/,$dp[3];
	my @dp2=split /=/,$dp[4];
	my @dp3=split /,/,$dp2[1];
	$genotype{0}=$line[3];
	my $i=1;
	foreach (@allel){
		$genotype{$i}=$_;
		$i++;
	}
	my$GT=$b[0];
	if(($GT!~/0\/0/ && $GT!~/\.\/\./) && (($dp1[1]<$HOMOZYGOUS) or ($b[-1]<$GQ) or ($line[5]<20))){
		print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t_\/_\n";
	}elsif($GT=~/0\/0/ && $line[5]>=20){
		print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$GT\n";
	}elsif($GT=~/\.\/\./){
		print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t_/_\n";
	}elsif($GT=~/1\/1/ && (($dp3[1]<$HOMOZYGOUS) or ($b[-1]<$GQ) or $line[5]<20)){
		print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t_/_\n";
	}elsif($GT=~/0\/1/ && (($dp3[0]<$HETEROZYGOUS) or ($dp3[1]<$HETEROZYGOUS) or ($b[-1]<$GQ) or $line[5]<20)){
		print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t_/_\n";
	}elsif($GT!~/0\/0/){
		$b[0]=~s:(\d).*(\d):$genotype{$1}/$genotype{$2}:g;
		print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$b[0]\n";
	}
		undef %genotype;
}
close IN;
close OU;
my %hash0;
my %tmp;
my %tmp1;
my %hash;
my @parent_id=glob("*_$inputfile\.parentgt");
my @Parent_Name=split/\_/,$parent_id[0];
my @progeny_id=glob("*_$inputfile\.gt");
open(IN,"01-$inputfile.vcf.gt");
while(<IN>){
	chop $_;
	if($_=~/#/){
		next;
	}else{
		my @fd = split/\s+/,$_;
		my $len1=length$fd[2];
		my $len2=($fd[1]+$len1-1);
		my $id="$fd[0]:$fd[1]";
		$hash0{$id} = "_/_";
	}
}
close IN;
open IN, "$parent_id[0]";
while(<IN>){
	chop $_;
	if($_=~/#/){
		next;
	}else{
		my @fd =split/\s+/,$_;
		my @fd1 = split/\-/,$fd[0];
		$hash0{$fd1[0]} = "_/_";
	}
}
close IN;
for(my $i=0;$i<$progeny_number;$i++){
	open IN,"$progeny_id[$i]";
	while(<IN>){
		chop $_;
		if($_=~/#/){
			next;
		}else{
			my @fd = split/\s+/,$_;
			my @fd1 = split/\-/,$fd[0];
			$hash0{$fd1[0]} = "_/_";
		}
	}
}
close IN;
open (IN,"<01-$inputfile.vcf.gt");
%hash=%hash0;
while(<IN>){
	chop$_;
	my @fd = split/\s+/,$_;
	my $len1=length$fd[2];
	my $len2=($fd[1]+$len1-1);
	my $id="$fd[0]:$fd[1]";
	$hash{$id} ="$fd[2]\t$fd[3]\t$fd[4]";
}
close(IN);
open IN, "$parent_id[0]";
%tmp=%hash0;
while(<IN>){
	chop$_;
	my @fd = split/\s+/,$_;
	my @fd1=split/\-/,$fd[0];
	$tmp{$fd1[0]} ="$fd[1]\t$fd[2]\t$fd[3]";
}
foreach my $key(keys %hash){
	$hash{$key}="$hash{$key}\t$tmp{$key}";
}
for(my $j=0;$j<$progeny_number;$j++){
	open IN,"$progeny_id[$j]";
	%tmp1 = %hash0;
	while(<IN>){
		chop$_;
		if($_=~/#/){
			next;
		}else{
			my @fd = split/\s+/,$_;
			my @fd1=split/\-/,$fd[0];
			$tmp1{$fd1[0]} = $fd[-1];
		}
	}
	close (IN);
	foreach my $key(keys %hash){
	$hash{$key}="$hash{$key}\t$tmp1{$key}";
	}
}
open OU,">total-$inputfile.loc";
foreach my $key(sort {(($a=~/(\d+)/)[0] <=> ($b=~/(\d+)/)[0]) or ((split/\:/,$a)[1] <=> (split/\:/,$b)[1])} (keys%hash)){
	print OU "$key\t$hash{$key}\n";
}
close OU;
undef%hash0;
undef%hash;
undef%tmp1;
open IN, "total-$inputfile.loc";
open OU,">total-$inputfile.loc1";
while(<IN>){
	chomp;
	my @line=split/\s+/,$_;
	if($line[1]!~/_\/_/){
		print OU "$_\n";
	}else{
		next;
	}
}
close IN;
close OU;
unlink "total-$inputfile.loc";
rename "total-$inputfile.loc1", "total-$inputfile.loc";
unlink "01-$inputfile.vcf";
unlink "01-$inputfile.vcf.gt";
my $t="@progeny_id";
$t=~s/ /\t/g;
$t=~s/\_progeny_$inputfile\.gt//g;
open SEIN, "total-$inputfile.loc";
open SEOU, ">$inputfile\_heterozygote.txt";
print SEOU "Position\t$inputfile\t$Parent_Name[0]\t$t\n";
while(<SEIN>){
	chomp;
	my @segregation=split/\s+/,$_;
	if(($segregation[1]=~/_\/_/) or ($segregation[3]=~/_\/_/) or ($segregation[4]=~/_\/_/) or ($segregation[6]=~/_\/_/)){
		next;
	}else{
		my $se1=shift@segregation;
		my $se2=shift@segregation;
		my $se3=shift@segregation;
		my $se4=shift@segregation;
		my $se5=shift@segregation;
		my $se6=shift@segregation;
		my $se7=shift@segregation;
		my @gt4=split/\//,$se4;
		my @gt7=split/\//,$se7;
		my $se_gt="@segregation";
		$se_gt=~s/\s/\t/g;
		$se_gt=~s/0\/0/$se2\/$se2/g;
		if($se4 eq $se7){
			print SEOU "$se1\t$se4(ab)\t$se7(ab)\t$se_gt\n";
		}elsif(($se4 ne $se7) && ($gt7[0] eq $gt7[1]) && ($se5 eq $se2) && ($gt7[0] eq "0")){
			print SEOU "$se1\t$se4(ab)\t$se2\/$se2(aa)\t$se_gt\n";
		}elsif(($se4 ne $se7) && ($gt7[0] eq $gt7[1]) && ($gt7[0] eq $gt4[1])){
			print SEOU "$se1\t$se4(ab)\t$se7(aa)\t$se_gt\n";
		}elsif(($se4 ne $se7) && ($gt7[0] eq $gt7[1]) && ($gt7[0] ne $gt4[0]) && ($gt7[0] ne $gt4[1])){
			print SEOU "$se1\t$se4(ab)\t$se7(cc)\t$se_gt\n";
		}elsif(($se4 ne $se7) && ($gt7[0] ne $gt7[1]) && ($gt7[0] eq $gt4[0])){
			print SEOU "$se1\t$se4(ab)\t$se7(ac)\t$se_gt\n";
		}elsif(($se4 ne $se7) && ($gt7[0] ne $gt7[1]) && ($gt7[0] eq $gt4[1])){
			print SEOU "$se1\t$gt4[1]/$gt4[0](ab)\t$se7(ac)\t$se_gt\n";
		}elsif(($se4 ne $se7) && ($gt7[0] ne $gt7[1]) && ($gt7[1] eq $gt4[0])){
			print SEOU "$se1\t$se4(ab)\t$gt7[1]/$gt7[0](ac)\t$se_gt\n";
		}elsif(($se4 ne $se7) && ($gt7[0] ne $gt7[1]) && ($gt7[1] eq $gt4[1])){
			print SEOU "$se1\t$gt4[1]/$gt4[0](ab)\t$gt7[1]/$gt7[0](ac)\t$se_gt\n";
		}elsif(($se4 ne $se7) && ($gt7[0] ne $gt7[1]) && ($gt7[0] ne $gt4[0]) && ($gt7[0] ne $gt4[1]) && ($gt7[1] ne $gt4[0]) && ($gt7[1] ne $gt4[1])){
			print SEOU "$se1\t$se4(ab)\t$se7(cd)\t$se_gt\n";
		}
	}
}
close SEIN;
close SEOU;
unlink "total-$inputfile.loc";
open ABIN, "$inputfile\_heterozygote.txt";
open ABXAA, ">$inputfile\_abxaa_gt.txt";
open ABXAB, ">abxab_gt.txt";
open ABXCC, ">$inputfile\_abxcc_gt.txt";
open ABXAC, ">abxac_gt.txt";
open ABXCD, ">abxcd_gt.txt";
while(<ABIN>){
        chomp;
	my @line=split/\s+/,$_;
        if($_=~/Position/){
                next;
	}elsif($_=~/aa/){
		my $id1=shift@line;
                my $id2=shift@line;
                my $id3=shift@line;
                my @gt1=split/\(/,$id2;
                my @gt2=split/\(/,$id3;
                my @abxaa_gt1=split/\//,$gt1[0];
                my @abxaa_gt2=split/\//,$gt2[0];
                my $nn="@line";
                my $n1=($nn=~s/\b$abxaa_gt1[0]\/$abxaa_gt1[1]\b/\b$abxaa_gt1[0]\/$abxaa_gt1[1]\b/g);
                my $n2=($nn=~s/\b$abxaa_gt2[0]\/$abxaa_gt2[1]\b/\b$abxaa_gt2[0]\/$abxaa_gt2[1]\b/g);
                my $n3=($nn=~s/_\/_/_\/_/g);
		my $n4=($nn=~s/\b$abxaa_gt1[1]\/$abxaa_gt1[0]\b/\b$abxaa_gt1[1]\/$abxaa_gt1[0]\b/g);
		my $n5=($n1+$n4);
                if($n5+$n2+$n3==$progeny_number && $n5!=0 && $n2!=0 && $n3==0){
                        print ABXAA "$id1\t$n5\t$n2\n";
                }elsif($n5+$n2+$n3==$progeny_number && $n5!=0 && $n2!=0 && $n3!=0 && ($n3/$progeny_number)<=($MISPCT/100)){
                        print ABXAA "$id1\t$n5\t$n2\n";
                }
	}elsif($line[2]=~/ab/){
		my $id1=shift@line;
                my $id2=shift@line;
                my $id3=shift@line;
                my @gt1=split/\(/,$id2;
                my @gt2=split/\(/,$id3;
                my @abxab_gt1=split/\//,$gt1[0];
                my @abxab_gt2=split/\//,$gt2[0];
                my $nn="@line";
                my $n1=($nn=~s/\b$abxab_gt1[0]\/$abxab_gt1[1]\b/\b$abxab_gt1[0]\/$abxab_gt1[1]\b/g);
                my $n2=($nn=~s/\b$abxab_gt1[0]\/$abxab_gt1[0]\b/\b$abxab_gt1[0]\/$abxab_gt1[0]\b/g);
		my $n3=($nn=~s/\b$abxab_gt1[1]\/$abxab_gt1[1]\b/\b$abxab_gt1[1]\/$abxab_gt1[1]\b/g);
                my $n4=($nn=~s/_\/_/_\/_/g);
		my $n5=($nn=~s/\b$abxab_gt1[1]\/$abxab_gt1[0]\b/\b$abxab_gt1[1]\/$abxab_gt1[0]\b/g);
		my $n6=($n1+$n5);
                if($n6+$n2+$n3+$n4==$progeny_number && $n6!=0 && $n2!=0 && $n3!=0 && $n4==0){
                        print ABXAB "$id1\t$n6\t$n2\t$n3\n";
                }elsif($n6+$n2+$n3+$n4==$progeny_number && $n6!=0 && $n2!=0 && $n3!=0 && $n4!=0 && ($n4/$progeny_number)<=($MISPCT/100)){
                        print ABXAB "$id1\t$n6\t$n2\t$n3\n";
                }
        }elsif($_=~/cc/){
                my $id1=shift@line;
                my $id2=shift@line;
                my $id3=shift@line;
                my @gt1=split/\(/,$id2;
                my @gt2=split/\(/,$id3;
		my @abxcc_gt1=split/\//,$gt1[0];
		my @abxcc_gt2=split/\//,$gt2[0];
                my $nn="@line";
                my $n1=($nn=~s/\b$abxcc_gt1[0]\/$abxcc_gt2[0]\b/\b$abxcc_gt1[0]\/$abxcc_gt2[0]\b/g);
                my $n2=($nn=~s/\b$abxcc_gt1[1]\/$abxcc_gt2[0]\b/\b$abxcc_gt1[1]\/$abxcc_gt2[0]\b/g);
		my $n3=($nn=~s/\b$abxcc_gt2[0]\/$abxcc_gt1[0]\b/\b$abxcc_gt2[0]\/$abxcc_gt1[0]\b/g);
		my $n4=($nn=~s/\b$abxcc_gt2[0]\/$abxcc_gt1[1]\b/\b$abxcc_gt2[0]\/$abxcc_gt1[1]\b/g);
		my $n5=($n1+$n3);
		my $n6=($n2+$n4);
                my $n7=($nn=~s/_\/_/_\/_/g);
                if($n5+$n6+$n7==$progeny_number && $n5!=0 && $n6!=0 && $n7==0){
                        print ABXCC "$id1\t$n5\t$n6\n";
                }elsif($n5+$n6+$n7==$progeny_number && $n5!=0 && $n6!=0 && $n7!=0 && ($n7/$progeny_number)<=($MISPCT/100)){
                        print ABXCC "$id1\t$n5\t$n6\n";
                }
        }elsif($_=~/ac/){
		my $id1=shift@line;
		my $id2=shift@line;
		my $id3=shift@line;
		my @gt1=split/\(/,$id2;
		my @gt2=split/\(/,$id3;
		my @abxac_gt1=split/\//,$gt1[0];
		my @abxac_gt2=split/\//,$gt2[0];
		my $nn="@line";
		$nn=~s/\b$abxac_gt1[1]\/$abxac_gt1[0]\b/\b$abxac_gt1[0]\/$abxac_gt1[1]\b/g;
		$nn=~s/\b$abxac_gt2[1]\/$abxac_gt2[0]\b/\b$abxac_gt2[0]\/$abxac_gt2[1]\b/g;
		$nn=~s/\b$abxac_gt2[1]\/$abxac_gt1[1]\b/\b$abxac_gt1[1]\/$abxac_gt2[1]\b/g;
		$nn=~s/\b$abxac_gt2[0]\/$abxac_gt1[0]\b/\b$abxac_gt1[0]\/$abxac_gt2[0]\b/g;
		my $n1=($nn=~s/\b$abxac_gt1[0]\/$abxac_gt1[1]\b/\b$abxac_gt1[0]\/$abxac_gt1[1]\b/g);
		my $n2=($nn=~s/\b$abxac_gt1[0]\/$abxac_gt1[0]\b/\b$abxac_gt1[0]\/$abxac_gt1[0]\b/g);
		my $n3=($nn=~s/\b$abxac_gt2[0]\/$abxac_gt2[1]\b/\b$abxac_gt2[0]\/$abxac_gt2[1]\b/g);
		my $n4=($nn=~s/\b$abxac_gt1[1]\/$abxac_gt2[1]\b/\b$abxac_gt1[1]\/$abxac_gt2[1]\b/g);
		my $n5=($nn=~s/_\/_/_\/_/g);
		if($n1+$n2+$n3+$n4+$n5==$progeny_number && $n1!=0 && $n2!=0 && $n3!=0 && $n4!=0 && $n5==0){
			print ABXAC "$id1\t$n1\t$n2\t$n3\t$n4\n";
		}elsif($n1+$n2+$n3+$n4+$n5==$progeny_number && $n1!=0 && $n2!=0 && $n3!=0 && $n4!=0 && $n5!=0 && ($n5/$progeny_number)<=($MISPCT/100)){
			print ABXAC "$id1\t$n1\t$n2\t$n3\t$n4\n";
		}
	}elsif($_=~/cd/){
                my $id1=shift@line;
                my $id2=shift@line;
                my $id3=shift@line;
                my @gt1=split/\(/,$id2;
                my @gt2=split/\(/,$id3;
                my @abxcd_gt1=split/\//,$gt1[0];
                my @abxcd_gt2=split/\//,$gt2[0];
                my $nn="@line";
		$nn=~s/\b$abxcd_gt2[0]\/$abxcd_gt1[0]\b/\b$abxcd_gt1[0]\/$abxcd_gt2[0]\b/g;
		$nn=~s/\b$abxcd_gt2[1]\/$abxcd_gt1[0]\b/\b$abxcd_gt1[0]\/$abxcd_gt2[1]\b/g;
		$nn=~s/\b$abxcd_gt2[0]\/$abxcd_gt1[1]\b/\b$abxcd_gt1[1]\/$abxcd_gt2[0]\b/g;
		$nn=~s/\b$abxcd_gt2[1]\/$abxcd_gt1[1]\b/\b$abxcd_gt1[1]\/$abxcd_gt2[1]\b/g;
		my $n1=($nn=~s/\b$abxcd_gt1[0]\/$abxcd_gt2[0]\b/\b$abxcd_gt1[0]\/$abxcd_gt2[0]\b/g);
		my $n2=($nn=~s/\b$abxcd_gt1[0]\/$abxcd_gt2[1]\b/\b$abxcd_gt1[0]\/$abxcd_gt2[1]\b/g);
		my $n3=($nn=~s/\b$abxcd_gt1[1]\/$abxcd_gt2[0]\b/\b$abxcd_gt1[1]\/$abxcd_gt2[0]\b/g);
		my $n4=($nn=~s/\b$abxcd_gt1[1]\/$abxcd_gt2[1]\b/\b$abxcd_gt1[1]\/$abxcd_gt2[1]\b/g);
		my $n5=($nn=~s/_\/_/_\/_/g);
		if($n1+$n2+$n3+$n4+$n5==$progeny_number && $n1!=0 && $n2!=0 && $n3!=0 && $n4!=0 && $n5==0){
			print ABXCD "$id1\t$n1\t$n2\t$n3\t$n4\n";
		}elsif($n1+$n2+$n3+$n4+$n5==$progeny_number && $n1!=0 && $n2!=0 && $n3!=0 && $n4!=0 && $n5!=0 && ($n5/$progeny_number)<=($MISPCT/100)){
			print ABXCD "$id1\t$n1\t$n2\t$n3\t$n4\n";
		}
	}
}
close ABIN;
close ABXAA;
close ABXAB;
close ABXCC;
close ABXAC;
close ABXCD;
chdir '..';
my %ab_pchisq_segregation;
my @abxaa_line1;
my @abxcc_line1;
my @abxab_line1;
my @abxac_line1;
my @abxcd_line1;
my $abxaa_line=(`wc -l $directory/$inputfile\_abxaa_gt.txt`);
my $abxcc_line=(`wc -l $directory/$inputfile\_abxcc_gt.txt`);
my $abxab_line=(`wc -l $directory/abxab_gt.txt`);
my $abxac_line=(`wc -l $directory/abxac_gt.txt`);
my $abxcd_line=(`wc -l $directory/abxcd_gt.txt`);
if($population eq 'CP'){
	@abxaa_line1=split/\s+/,$abxaa_line;
	@abxcc_line1=split/\s+/,$abxcc_line;
	@abxab_line1=split/\s+/,$abxab_line;
	@abxac_line1=split/\s+/,$abxac_line;
	@abxcd_line1=split/\s+/,$abxcd_line;
}
if($population eq 'F2'){
        @abxaa_line1=split/\s+/,$abxaa_line;
        @abxcc_line1=split/\s+/,$abxcc_line;
        @abxab_line1=split/\s+/,$abxab_line;
        @abxac_line1=split/\s+/,$abxac_line;
        @abxcd_line1=split/\s+/,$abxcd_line;
	$abxaa_line1[0]=0;
	$abxcc_line1[0]=0;
	$abxac_line1[0]=0;
	$abxcd_line1[0]=0;
}
if($population eq 'BC1' or $population eq 'BC2'){
        @abxaa_line1=split/\s+/,$abxaa_line;
        @abxcc_line1=split/\s+/,$abxcc_line;
        @abxab_line1=split/\s+/,$abxab_line;
        @abxac_line1=split/\s+/,$abxac_line;
        @abxcd_line1=split/\s+/,$abxcd_line;
        $abxab_line1[0]=0;
        $abxcc_line1[0]=0;
        $abxac_line1[0]=0;
        $abxcd_line1[0]=0;
}

if($abxcc_line1[0]!=0){
	open CCPVIN, "$directory/$inputfile\_abxcc_gt.txt";
	open CCPVOU, ">$directory/$inputfile-abxcc-pchisq.segregation_ratio";
	while(<CCPVIN>){
		chomp;
		my @ccpv=split/\s+/,$_;
		my $ccmean=(($ccpv[1]+$ccpv[2])/2);
		my $ccx2=($ccpv[1]-$ccmean)**2/$ccmean + ($ccpv[2]-$ccmean)**2/$ccmean;
		my $ccpvalue=Statistics::Distributions::chisqrprob(1,$ccx2);
		if($ccpvalue>$PVALUE){
			print CCPVOU "$ccpv[0]\t$ccpv[1]\t$ccpv[2]\t$ccmean\t$ccpvalue\n";
		}else{
			next;
		}
	}
	close CCPVIN;
	close CCPVOU;
	open ABXCCSE,"$directory/$inputfile-abxcc-pchisq.segregation_ratio";
	while(<ABXCCSE>){
		chomp;
		my @abxcc_segregation_id=split/\s+/,$_;
		$ab_pchisq_segregation{$abxcc_segregation_id[0]}="$abxcc_segregation_id[1]\t$abxcc_segregation_id[2]\t$abxcc_segregation_id[4]";
	}
	open ABXCCIN,"$directory/$inputfile\_heterozygote.txt";
	open ABXCCOUT,">$directory/$inputfile\_abxcc.txt";
	print ABXCCOUT "Position\t$inputfile\t$Parent_Name[0]\t$t\tthe_number_of_ac\tthe_number_of_bc\tthe_number_of_missing\tPVALUE\n";
	while(<ABXCCIN>){
		chomp;
		my @abxcc_p=split/\s+/,$_;
		if(exists $ab_pchisq_segregation{$abxcc_p[0]}){
			my @abcc1=split/\(/,$abxcc_p[1];
			my @abcc2=split/\(/,$abxcc_p[2];
                        my @abxcc_gt1=split/\//,$abcc1[0];
                        my @abxcc_gt2=split/\//,$abcc2[0];
                        $_=~s/\(ab\)//g;
                        $_=~s/\(cc\)//g;
                        $_=~s/\b$abcc1[0]\b/ab/g;
                        $_=~s/\b$abcc2[0]\b/cc/g;
                        $_=~s/\b$abxcc_gt2[0]\/$abxcc_gt1[0]\b/ac/g;
                        $_=~s/\b$abxcc_gt2[0]\/$abxcc_gt1[1]\b/bc/g;
			$_=~s/\b$abxcc_gt1[0]\/$abxcc_gt2[0]\b/ac/g;
			$_=~s/\b$abxcc_gt1[1]\/$abxcc_gt2[0]\b/bc/g;
			$_=~s/_\/_/--/g;
			my @abxcc_missing_num=split/\s+/,$ab_pchisq_segregation{$abxcc_p[0]};
			my $abccmiss=($progeny_number-($abxcc_missing_num[0]+$abxcc_missing_num[1]));
			print ABXCCOUT "$_\tac:$abxcc_missing_num[0]\tbc:$abxcc_missing_num[1]\t--:$abccmiss\t$abxcc_missing_num[2]\n";
		}else{
			next;
		}
	}
	undef %ab_pchisq_segregation;
	unlink "$directory/$inputfile\_abxcc_gt.txt","$directory/$inputfile-abxcc-pchisq.segregation_ratio";
}elsif($abxcc_line1[0]==0){
	unlink "$directory/$inputfile\_abxcc_gt.txt";
}
if($abxac_line1[0]==0){
        unlink "$directory/abxac_gt.txt";
}elsif($abxac_line1[0]!=0){
        open ACPVIN, "$directory/abxac_gt.txt";
        open ACPVOU, ">$directory/abxac-pchisq.segregation_ratio";
        while(<ACPVIN>){
                chomp;
                my @acpv=split/\s+/,$_;
                my $acmean=(($acpv[1]+$acpv[2]+$acpv[3]+$acpv[4])/4);
                my $acx2=($acpv[1]-$acmean)**2/$acmean + ($acpv[2]-$acmean)**2/$acmean + ($acpv[3]-$acmean)**2/$acmean + ($acpv[4]-$acmean)**2/$acmean;
                my $acpvalue=Statistics::Distributions::chisqrprob(3,$acx2);
		if($acpvalue>$PVALUE){
	                print ACPVOU "$acpv[0]\t$acpv[1]\t$acpv[2]\t$acpv[3]\t$acpv[4]\t$acmean\t$acpvalue\n";
		}else{
			next;
		}
        }
	close ACPVIN;
	close ACPVOU;
	open ABXACSE,"$directory/abxac-pchisq.segregation_ratio";
	while(<ABXACSE>){
		chomp;
		my @abxac_segregation_id=split/\s+/,$_;
		$ab_pchisq_segregation{$abxac_segregation_id[0]}="$abxac_segregation_id[1]\t$abxac_segregation_id[2]\t$abxac_segregation_id[3]\t$abxac_segregation_id[4]\t$abxac_segregation_id[6]";
	}
	open ABXACIN,"$directory/$inputfile\_heterozygote.txt";
	open ABXACOUT,">$directory/abxac.txt";
	print ABXACOUT "Position\t$inputfile\t$Parent_Name[0]\t$t\tthe_number_of_ab\tthe_number_of_aa\tthe_number_of_ac\tthe_number_of_bc\tthe_number_of_missing\tPVALUE\n";
	while(<ABXACIN>){
		chomp;
		my @abxac_p=split/\s+/,$_;
		if(exists $ab_pchisq_segregation{$abxac_p[0]}){
                        my @abac1=split/\(/,$abxac_p[1];
                        my @abac2=split/\(/,$abxac_p[2];			
			my @abxac_gt1=split/\//,$abac1[0];
			my @abxac_gt2=split/\//,$abac2[0];
	                $_=~s/\b$abxac_gt1[1]\/$abxac_gt1[0]\b/ab/g;
                	$_=~s/\b$abxac_gt2[1]\/$abxac_gt2[0]\b/ac/g;
			$_=~s/\b$abxac_gt1[0]\/$abxac_gt2[0]\b/aa/g;
			$_=~s/\b$abxac_gt2[1]\/$abxac_gt1[1]\b/bc/g;
			$_=~s/\b$abxac_gt1[1]\/$abxac_gt2[1]\b/bc/g;
			$_=~s/\b$abac1[0]\b/ab/g;
			$_=~s/\b$abac2[0]\b/ac/g;
			$_=~s/\(ab\)//g;
			$_=~s/\(ac\)//g;
			$_=~s/_\/_/--/g;
			my @abxac_missing_num=split/\s+/,$ab_pchisq_segregation{$abxac_p[0]};
			my $abacmiss=($progeny_number-($abxac_missing_num[0]+$abxac_missing_num[1]+$abxac_missing_num[2]+$abxac_missing_num[3]));
			print ABXACOUT "$_\tab:$abxac_missing_num[0]\taa:$abxac_missing_num[1]\tac:$abxac_missing_num[2]\tbc:$abxac_missing_num[3]\t--:$abacmiss\t$abxac_missing_num[4]\n";
		}else{
			next;
		}
	}
	undef %ab_pchisq_segregation;
	unlink "$directory/abxac_gt.txt","$directory/abxac-pchisq.segregation_ratio";
}
if($abxcd_line1[0]==0){
        unlink "$directory/abxcd_gt.txt";
}elsif($abxcd_line1[0]!=0){
        open CDPVIN, "$directory/abxcd_gt.txt";
        open CDPVOU, ">$directory/abxcd-pchisq.segregation_ratio";
        while(<CDPVIN>){
                chomp;
                my @cdpv=split/\s+/,$_;
                my $cdmean=(($cdpv[1]+$cdpv[2]+$cdpv[3]+$cdpv[4])/4);
                my $cdx2=($cdpv[1]-$cdmean)**2/$cdmean + ($cdpv[2]-$cdmean)**2/$cdmean + ($cdpv[3]-$cdmean)**2/$cdmean + ($cdpv[4]-$cdmean)**2/$cdmean;
                my $cdpvalue=Statistics::Distributions::chisqrprob(3,$cdx2);
		if($cdpvalue>$PVALUE){
	                print CDPVOU "$cdpv[0]\t$cdpv[1]\t$cdpv[2]\t$cdpv[3]\t$cdpv[4]\t$cdmean\t$cdpvalue\n";
		}else{
			next;
		}
        }
	close CDPVIN;
	close CDPVOU;
	open ABXCDSE,"$directory/abxcd-pchisq.segregation_ratio";
        while(<ABXCDSE>){
                chomp;
                my @abxcd_segregation_id=split/\s+/,$_;
                $ab_pchisq_segregation{$abxcd_segregation_id[0]}="$abxcd_segregation_id[1]\t$abxcd_segregation_id[2]\t$abxcd_segregation_id[3]\t$abxcd_segregation_id[4]\t$abxcd_segregation_id[6]";
        }
        open ABXCDIN,"$directory/$inputfile\_heterozygote.txt";
        open ABXCDOUT,">$directory/abxcd.txt";
	print ABXCDOUT "Position\t$inputfile\t$Parent_Name[0]\t$t\tthe_number_of_ac\tthe_number_of_ad\tthe_number_of_bc\tthe_number_of_bd\tthe_number_of_missing\tPVALUE\n";
        while(<ABXCDIN>){
                chomp;
                my @abxcd_p=split/\s+/,$_;
                if(exists $ab_pchisq_segregation{$abxcd_p[0]}){
                        my @abcd1=split/\(/,$abxcd_p[1];
                        my @abcd2=split/\(/,$abxcd_p[2];
                        my @abxcd_gt1=split/\//,$abcd1[0];
                        my @abxcd_gt2=split/\//,$abcd2[0];
			$_=~s/\(ab\)//g;
			$_=~s/\(cd\)//g;
			$_=~s/$abcd1[0]/ab/g;
			$_=~s/$abcd2[0]/cd/g;
			$_=~s/\b$abxcd_gt2[0]\/$abxcd_gt1[0]\b/ac/g;
	                $_=~s/\b$abxcd_gt2[1]\/$abxcd_gt1[0]\b/ad/g;
        	        $_=~s/\b$abxcd_gt2[0]\/$abxcd_gt1[1]\b/bc/g;
                	$_=~s/\b$abxcd_gt2[1]\/$abxcd_gt1[1]\b/bd/g;
			$_=~s/\b$abxcd_gt1[0]\/$abxcd_gt2[0]\b/ac/g;
                        $_=~s/\b$abxcd_gt1[0]\/$abxcd_gt2[1]\b/ad/g;
                        $_=~s/\b$abxcd_gt1[1]\/$abxcd_gt2[0]\b/bc/g;
                        $_=~s/\b$abxcd_gt1[1]\/$abxcd_gt2[1]\b/bd/g;
			$_=~s/_\/_/--/g;
                        my @abxcd_missing_num=split/\s+/,$ab_pchisq_segregation{$abxcd_p[0]};
                        my $abcdmiss=($progeny_number-($abxcd_missing_num[0]+$abxcd_missing_num[1]+$abxcd_missing_num[2]+$abxcd_missing_num[3]));
                        print ABXCDOUT "$_\tac:$abxcd_missing_num[0]\tad:$abxcd_missing_num[1]\tbc:$abxcd_missing_num[2]\tbd:$abxcd_missing_num[3]\t--:$abcdmiss\t$abxcd_missing_num[4]\n";
                }else{
                        next;
                }
        }
	undef %ab_pchisq_segregation;
	unlink "$directory/abxcd_gt.txt","$directory/abxcd-pchisq.segregation_ratio";
}
if($abxaa_line1[0]==0){
	unlink "$directory/$inputfile\_abxaa_gt.txt";
}elsif($abxaa_line1[0]!=0){
	open AAPVIN, "$directory/$inputfile\_abxaa_gt.txt";
	open AAPVOU, ">$directory/$inputfile-abxaa-pchisq.segregation_ratio";
	while(<AAPVIN>){
		chomp;
		my @aapv=split/\s+/,$_;
		my $aamean=(($aapv[1]+$aapv[2])/2);
		my $aax2=($aapv[1]-$aamean)**2/$aamean + ($aapv[2]-$aamean)**2/$aamean;
		my $aapvalue=Statistics::Distributions::chisqrprob(1,$aax2);
		if($aapvalue>$PVALUE){
			print AAPVOU "$aapv[0]\t$aapv[1]\t$aapv[2]\t$aamean\t$aapvalue\n";
		}else{
			next;
		}
	}
	close AAPVIN;
	close AAPVOU;
	open ABXAASE,"$directory/$inputfile-abxaa-pchisq.segregation_ratio";
	while(<ABXAASE>){
		chomp;
		my @abxaa_segregation_id=split/\s+/,$_;
		$ab_pchisq_segregation{$abxaa_segregation_id[0]}="$abxaa_segregation_id[1]\t$abxaa_segregation_id[2]\t$abxaa_segregation_id[4]";
	}
	open ABXAAIN,"$directory/$inputfile\_heterozygote.txt";
	open ABXAAOUT,">$directory/$inputfile\_abxaa.txt";
	print ABXAAOUT "Position\t$inputfile\t$Parent_Name[0]\t$t\tthe_number_of_ab\tthe_number_of_aa\tthe_number_of_missing\tPVALUE\n";
		while(<ABXAAIN>){
		chomp;
		my @abxaa_p=split/\s+/,$_;
		if(exists $ab_pchisq_segregation{$abxaa_p[0]}){
			my @abaa1=split/\(/,$abxaa_p[1];
			my @abaa2=split/\(/,$abxaa_p[2];
			my @abxaa_gt1=split/\//,$abaa1[0];
			my @abxaa_gt2=split/\//,$abaa2[0];
			$_=~s/\(ab\)//g;
			$_=~s/\(aa\)//g;
			$_=~s/\b$abxaa_gt1[1]\/$abxaa_gt1[0]\b/ab/g;
			$_=~s/\b$abxaa_gt1[0]\/$abxaa_gt1[1]\b/ab/g;
			$_=~s/\b$abxaa_gt2[0]\/$abxaa_gt2[0]\b/aa/g;
			$_=~s/_\/_/--/g;
			my @abxaa_missing_num=split/\s+/,$ab_pchisq_segregation{$abxaa_p[0]};
			my $abaamiss=($progeny_number-($abxaa_missing_num[0]+$abxaa_missing_num[1]));
			print ABXAAOUT "$_\tab:$abxaa_missing_num[0]\taa:$abxaa_missing_num[1]\t--:$abaamiss\t$abxaa_missing_num[2]\n";
		}else{
			next;
		}
	}
	undef %ab_pchisq_segregation;
	unlink "$directory/$inputfile\_abxaa_gt.txt","$directory/$inputfile\-abxaa-pchisq.segregation_ratio";
}
if($abxab_line1[0]==0){
	unlink "$directory/abxab_gt.txt";
}elsif($abxab_line1[0]!=0){
        open ABPVIN, "$directory/abxab_gt.txt";
        open ABPVOU, ">$directory/abxab-pchisq.segregation_ratio";
        while(<ABPVIN>){
                chomp;
                my @abpv=split/\s+/,$_;
                my $abmean=(($abpv[1]+$abpv[2]+$abpv[3])/4);
                my $abx2=($abpv[1]-2*$abmean)**2/(2*$abmean) + ($abpv[2]-$abmean)**2/$abmean + ($abpv[3]-$abmean)**2/$abmean;
                my $abpvalue=Statistics::Distributions::chisqrprob(2,$abx2);
                if($abpvalue>$PVALUE){
                        print ABPVOU "$abpv[0]\t$abpv[1]\t$abpv[2]\t$abpv[3]\t$abmean\t$abpvalue\n";
                }else{
                        next;
                }
        }
        close ABPVIN;
        close ABPVOU;
	open ABXABSE,"$directory/abxab-pchisq.segregation_ratio";
	while(<ABXABSE>){
		chomp;
		my @abxab_segregation_id=split/\s+/,$_;
		$ab_pchisq_segregation{$abxab_segregation_id[0]}="$abxab_segregation_id[1]\t$abxab_segregation_id[2]\t$abxab_segregation_id[3]\t$abxab_segregation_id[5]";
	}
	open ABXABIN,"$directory/$inputfile\_heterozygote.txt";
	open ABXABOUT,">$directory/abxab.txt";
	print ABXABOUT "Position\t$inputfile\t$Parent_Name[0]\t$t\tthe_number_of_ab\tthe_number_of_aa\tthe_number_of_bb\tthe_number_of_missing\tPVALUE\n";
	while(<ABXABIN>){
		chomp;
		my @abxab_p=split/\s+/,$_;
		if(exists $ab_pchisq_segregation{$abxab_p[0]}){
		my @abab1=split/\(/,$abxab_p[1];
		my @abab2=split/\(/,$abxab_p[2];
		my @abxab_gt1=split/\//,$abab1[0];
		my @abxab_gt2=split/\//,$abab2[0];
		$_=~s/\(ab\)//g;
		$_=~s/\b$abxab_gt1[1]\/$abxab_gt1[0]\b/ab/g;
		$_=~s/\b$abxab_gt1[0]\/$abxab_gt1[1]\b/ab/g;
		$_=~s/\b$abxab_gt1[0]\/$abxab_gt1[0]\b/aa/g;
		$_=~s/\b$abxab_gt1[1]\/$abxab_gt1[1]\b/bb/g;
		$_=~s/_\/_/--/g;
		my @abxab_missing_num=split/\s+/,$ab_pchisq_segregation{$abxab_p[0]};
		my $ababmiss=($progeny_number-($abxab_missing_num[0]+$abxab_missing_num[1]+$abxab_missing_num[2]));
		print ABXABOUT "$_\tab:$abxab_missing_num[0]\taa:$abxab_missing_num[1]\tbb:$abxab_missing_num[2]\t--:$ababmiss\t$abxab_missing_num[3]\n";
		}else{
			next;
		}
	}
	undef %ab_pchisq_segregation;
	unlink "$directory/abxab_gt.txt","$directory/abxab-pchisq.segregation_ratio";
}
