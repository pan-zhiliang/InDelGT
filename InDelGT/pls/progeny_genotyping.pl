#!/usr/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;

sub prtHelp{
	print "\nThis procedure is used for indels genotyping of F1 generation populations\n";
	print "Contact: Chunfa Tong  <tongchf\@njfu.edu.cn>\n";
	print "Usage: perl genotyping.pl [Options] <*.vcf> [Options] <*.bcf>\n";
	print "\n";
	print "Options:\n";
	print "	-q  <int>	the minimum score of InDel genotypes\n";
	print "	-c  <str>	BCFtools fold\n";
	print "	-ho  <int>	the depth of homozygate\n";
	print "	-he  <int>	the depth of heterozygate\n";
	print "	-v  <str>	the catalog input file for one of the parent\n";
	print "	-b  <str>	the BCF input file of progeny\n";
	print "	-o  <str>	create a directory for storing output file of the genotyping results\n";
	print "	--help|h	help\n";

}

my $inputfile;
my $progeny_inputfile;
my $outputfile;
my $help;
my $i;
my $j;
my %hash;
my %hash1;
my %prs;
my $bcftools;
my $GQ;
my $HETEROZYGOUS;
my $HOMOZYGOUS;
my ($str1,$str2);
my @parentname;
my @progenyname;
my $keys;
if(@ARGV<7){
	prtHelp();
	exit;
}

GetOptions(
	"q:s"=>\$GQ,
	"c:s"=>\$bcftools,
	"ho:s"=>\$HOMOZYGOUS,
	"he:s"=>\$HETEROZYGOUS,
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
	print STDERR "\nError: the VCF file of the parent was not provided!\n\n";
	prtHelp();
	exit;
}

unless($progeny_inputfile){
	print STDERR "\nError: the BCF file of another parent or progeny was not provided!\n\n";
	prtHelp();
	exit;
}

unless($outputfile){
	print STDERR "\nError: please create a directory for storing output file!\n\n";
	prtHelp();
	exit;
}

unless($bcftools){
        print STDERR "\nError: BCFtools fold was not provided!\n\n";
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
my $cwd=getcwd;
$outputfile=$cwd;
mkdir $outputfile or die "Error: can't create directory '$outputfile' : $!" unless( -d $outputfile);
chdir $outputfile or die "Error: can't cd to directory '$outputfile' : $!";
my $progeny_inputfilebcf=$progeny_inputfile;
my $inputfilevcf=$inputfile;
$progeny_inputfilebcf=~s/\.bcf//g;
$inputfilevcf=~s/\.cls//g;
open IN, "$inputfile";
open OU, ">$progeny_inputfilebcf-$inputfilevcf.vcf.loci";
open OUT, ">$progeny_inputfilebcf-$inputfilevcf.catalog";
while(<IN>){
        chomp;
	my @line=split/\s+/,$_;
	my $len1=length($line[2]);
	my $len2=($line[1]+$len1-1);
	print OU "$line[0]\t$line[1]\t$len2\n";
	print OUT "$line[0]\t$line[1]\t$len2\t$line[2]\n";
}
close IN;
close OU;
close OUT;
if(!-e "$progeny_inputfile.csi"){
	`$bcftools/bcftools index $progeny_inputfile`;
	`$bcftools/bcftools call -m -f gq -T $progeny_inputfilebcf-$inputfilevcf.vcf.loci $progeny_inputfile | awk '\$1!~/^#/' > $progeny_inputfilebcf-$inputfilevcf.vcf`;
}if(-e "$progeny_inputfile.csi"){
	unlink "$progeny_inputfile.csi";
	`$bcftools/bcftools index $progeny_inputfile`;
	`$bcftools/bcftools call -m -f gq -T $progeny_inputfilebcf-$inputfilevcf.vcf.loci $progeny_inputfile | awk '\$1!~/^#/' > $progeny_inputfilebcf-$inputfilevcf.vcf`;
}
unlink "$progeny_inputfilebcf-$inputfilevcf.vcf.loci";
open IN, "$progeny_inputfilebcf-$inputfilevcf.catalog";
while(<IN>){
        chomp;
        my @line=split/\s+/,$_;
        my $l="$line[0]\t$line[1]";
        $hash1{$l}="$line[0]\t$line[1]\t$line[3]";
}
open IN, "$progeny_inputfilebcf-$inputfilevcf.vcf";
while(<IN>){
	chomp;
	my @line=split/\s+/,$_;
	my $l="$line[0]\_$line[1]";
	$hash{$l}=$_;
}
close IN;
open IN,"$progeny_inputfilebcf-$inputfilevcf.catalog";
open OU, ">$progeny_inputfilebcf-$inputfilevcf.loc.vcf";
while(<IN>){
	chomp;
	my @line1=split/\s+/,$_;
	my $k="$line1[0]\_$line1[1]";
	if(exists $hash{$k} && $hash{$k}=~/INDEL/){
		print OU "$hash{$k}\n";
	}elsif(exists $hash{$k} && $hash{$k}=~/0\/0/){
		print OU "$hash{$k}\n";
	}else{
		print OU "$line1[0]\t$line1[1]\t$line1[3]\t\.\t_/_\n";
	}
}
close IN;
close OU;
undef %hash;
open IN,"$progeny_inputfilebcf-$inputfilevcf.loc.vcf";
open OU,">$progeny_inputfile.vcf.gt";
while(<IN>){
	chomp;
	if($_=~/_\/_/){
		print OU "$_\n";
	}else{
		my @line=split/\s+/,$_;
	        my @dp=split /;/,$line[7];
	        my @cd=split /:/,$line[-1];
	        my @dp1=split /=/,$dp[3];
	        my @dp2=split /=/,$dp[4];
	        my @dp3=split /,/,$dp2[1];
		my @b=split/\:/,$line[-1];
		my @gn=split/\//,$b[0];
		my $k= "$line[0]\t$line[1]";
		my @allel=split/,/,$line[4];
		$hash{0}=$line[3];
                my $i=1;
		foreach (@allel){
			$hash{$i}=$_;
			$i++;
		}
		my $GT=$b[0];
		if(($GT!~/0\/0/ && $GT!~/\.\/\./) && (($dp1[1]<$HOMOZYGOUS) or ($cd[-1]<$GQ) or ($line[5]<20))){
                        print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t_\/_\n";
                }elsif(exists $hash1{$k} && $GT=~/0\/0/ && $line[5]>=20){
                        print OU "$hash1{$k}\t$line[4]\t$GT\n";
		}elsif($GT=~/0\/0/ && $line[5]<20){
			print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t_/_\n";
                }elsif($GT=~/\.\/\./){
                        print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t_/_\n";
                }elsif($GT!~/0\/0/ && ($gn[0]!=$gn[1]) && (($dp3[$gn[0]]<$HETEROZYGOUS) or ($dp3[$gn[1]]<$HETEROZYGOUS) or ($cd[-1]<$GQ) or ($line[5]<20))){
                        print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t_/_\n";
                }elsif($GT!~/0\/0/ && ($gn[0]==$gn[1]) && (($dp3[$gn[0]]<$HOMOZYGOUS) or ($cd[-1]<$GQ) or ($line[5]<20))){
                        print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t_/_\n";
                }elsif($GT!~/0\/0/){
                        $b[0]=~s:(\d).*(\d):$hash{$1}/$hash{$2}:g;
                        print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$b[0]\n";
                }
	}
}
close IN;
close OU;
undef %hash;
undef %hash1;
open IN, "$progeny_inputfile.vcf.gt";
open OU, ">$progeny_inputfile.vcf.gt00";
while(<IN>){
	chomp;
	my @line=split/\s+/,$_;
	if($line[-1]=~/0\/0/){
		my $length=length($line[2]);
		my $l1=shift@line;
		my $l2=shift@line;
		my $n="@line";
		$n=~s/\s/\t/g;
		my $length1=($l2+$length-1);
		print OU "$l1:$l2-$length1\t$n\n";
	}else{
		next;
	}
}
open IN,"$progeny_inputfilebcf-$inputfilevcf.vcf";
open OU,">$progeny_inputfile\-01";
while(<IN>){
	chomp;
	my @line=split/\s+/,$_;
	if($_=~/#/){
		next;
	}elsif($line[5]<20){
		print OU "$_\n";
	}elsif($line[4]!~/\./ or $line[-1]=~/\.\/\./){
		print OU "$_\n";
	}elsif($_=~/0\/0/ && $_!~/INDEL/){
		my @line1=split/\;/,$line[7];
		my @line2=split/\=/,$line1[0];
		if($line2[1]<$HOMOZYGOUS && $line2[1]>=0){
			print OU "$_\n";
		}else{
			next;
		}
	}elsif($_=~/0\/0/ && $_=~/INDEL/){
		print OU "$_\n";
	}else{
		next;
	}
}
close IN;
close OU;
open IN,"$progeny_inputfile\-01";
while(<IN>){
        chomp;
        my @line=split/\s+/,$_;
	$hash{"$line[0]\t$line[1]"}=$_;
}
close IN;
open IN,"$progeny_inputfile.vcf.gt00";
open OU,">$progeny_inputfile.vcf.gt01";
open TXT,">$progeny_inputfile.vcf.gt02";
while(<IN>){
	chomp;
	my @line=split/\s+/,$_;
	my @id=split/\:/,$line[0];
	my @id1=split/\-/,$id[1];
	my $c=($id1[1]-$id1[0]);
	for (my $i=0;$i<=$c;$i++){
		my $No=($id1[0]+$i);
		if(exists $hash{"$id[0]\t$No"}){
			print OU "$_\n";
		}
		if(exists $hash{"$id[0]\t$No"} && $No!=$id1[1]){
			print TXT "$_\n";
		}
	}
}
close IN;
close OU;
close TXT;
open IN,"$progeny_inputfile.vcf.gt00";
open OU,">$progeny_inputfile.vcf.gt03";
while(<IN>){
        chomp;
        my @line=split/\s+/,$_;
        my @id=split/\:/,$line[0];
        my @id1=split/\-/,$id[1];
	my $fd1="$id[0]\t$id1[1]";
	if(exists $hash{$fd1} && $hash{$fd1}=~/INDEL/){
		print OU "$_\n";
	}else{
		next;
	}
}
undef %hash;
close IN;
close OU;
my @files=("$progeny_inputfile.vcf.gt01","$progeny_inputfile.vcf.gt02","$progeny_inputfile.vcf.gt03");
for(my $i=0;$i<3;$i++){
	`less -S $files[$i] |uniq > $files[$i]1`;
undef %hash;
close IN;
close OU;
unlink "$files[$i]";
`mv $files[$i]1 $files[$i]`;
}
open IN,"$progeny_inputfile.vcf.gt01";
while(<IN>){
	chomp;
	my @line=split/\s+/,$_;
	$hash{$line[0]}=1;
}	
open IN,"$progeny_inputfile.vcf.gt00";
open OU,">$progeny_inputfile.vcf.gt001";
while(<IN>){
	chomp;
	my @line=split/\s+/,$_;
	if(exists $hash{$line[0]}){
		next;
	}else{
		print OU "$_\n";
	}
}
undef %hash;
open IN,"$progeny_inputfile.vcf.gt02";
while(<IN>){
        chomp;
        my @line=split/\s+/,$_;
        $hash{$line[0]}=1;
}
open IN,"$progeny_inputfile.vcf.gt03";
open OU,">$progeny_inputfile.vcf.gt031";
while(<IN>){
        chomp;
        my @line=split/\s+/,$_;
        if(exists $hash{$line[0]}){
                next;
        }else{
                print OU "$_\n";
        }
}
undef %hash;
`cat $progeny_inputfile.vcf.gt001 $progeny_inputfile.vcf.gt031 | sort > $progeny_inputfile.homo`;
unlink glob "$progeny_inputfile.vcf.gt0*";
unlink glob "$progeny_inputfile\-0*";
open IN, "$progeny_inputfile.homo";
while(<IN>){
        chomp;
        my @line=split/\s+/,$_;
        my @line1=split/\-/,$line[0];
        $hash{$line1[0]}=1;
}
open IN, "$progeny_inputfile.vcf.gt";
open OU, ">$progeny_inputfilebcf\.gt";
while(<IN>){
        chomp;
        my @line=split/\s+/,$_;
        my $length=length($line[2]);
        my $l1=shift@line;
        my $l2=shift@line;
        my $n="@line";
        $n=~s/\s/\t/g;
        my $length1=($l2+$length-1);
        my $m="$l1\:$l2";
        if(exists $hash{$m}){
                print OU "$l1:$l2-$length1\t$n\n";
        }elsif($n=~/0\/0/){
                $n=~s/0\/0/_\/_/g;
                print OU "$l1:$l2-$length1\t$n\n";
        }else{
                print OU "$l1:$l2-$length1\t$n\n";
        }
}
undef %hash;
close IN;
close OU;
unlink "$progeny_inputfile.homo","$progeny_inputfilebcf-$inputfilevcf.vcf","$progeny_inputfile.vcf.gt","$progeny_inputfilebcf-$inputfilevcf.vcf.loci","$progeny_inputfilebcf-$inputfilevcf.catalog","$progeny_inputfilebcf-$inputfilevcf.loc.vcf";

