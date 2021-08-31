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
	print "	-v  <file>	the VCF input file for one of the two parents\n";
	print "	-b  <file>	the BCF input file of progeny\n";
	print "	-o  <directory>	create a directory for storing output file of the genotyping results\n";
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
my $INDELGT;
my ($str1,$str2);
my @parentname;
my @progenyname;
my $keys;
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
open(PRS,"../parameters.ini") || die "\nError\! This program needs a parameter file,namely \"parameters.ini\".\n\n";
while(<PRS>){
        chomp;
        if($_=~/\:/ && ($_=~/BCFTOOLS_FOLD/ or $_=~/GQ/ or $_=~/HETEROZYGOUS_DEPTH/ or $_=~/HOMOZYGOUS_DEPTH/ or $_=~/INDELGT_FOLD/)){
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

if(exists($prs{"BCFTOOLS_FOLD"})){
        $bcftools = $prs{"BCFTOOLS_FOLD"};
}else{
        die "Error! Please check the line of BCFTOOLS_FOLD in \"parameters.ini\".\n";
}
if(exists($prs{"GQ"})){
	$GQ = $prs{"GQ"};
}else{
	die "Error! Please check the line of GQ in \"parameters.ini\".\n";
}
if(exists($prs{"HETEROZYGOUS_DEPTH"})){
        $HETEROZYGOUS = $prs{"HETEROZYGOUS_DEPTH"};
}else{
        die "Error! Please check the line of HETEROZYGOUS_DEPTH in \"parameters.ini\".\n";
}

if(exists($prs{"HOMOZYGOUS_DEPTH"})){
        $HOMOZYGOUS = $prs{"HOMOZYGOUS_DEPTH"};
}else{
        die "Error! Please check the line of HOMOZYGOUS_DEPTH in \"parameters.ini\".\n";
}
if(exists($prs{"INDELGT_FOLD"})){
        $INDELGT = $prs{"INDELGT_FOLD"};
}else{
        die "Error! Please check the line of INDELGT_FOLD in \"parameters.ini\".\n";
}

mkdir "$INDELGT/$outputfile" or die "Error: can't create directory '$outputfile' : $!" unless( -d "$INDELGT/$outputfile");
chdir "$INDELGT/$outputfile" or die "Error: can't cd to directory '$outputfile' : $!";
my $cwd=getcwd;
my $abs_inputfilepath="$cwd/$inputfile";
my $abs_progeny_inputfilepath="$cwd/$progeny_inputfile";
if($progeny_inputfile=~/\//){
	@progenyname=split/\//,$progeny_inputfile;
	$progenyname[-1]=~s/\.bcf//g;
}elsif($progeny_inputfile!~/\//){
	@progenyname=split/\//,$abs_progeny_inputfilepath;
	$progenyname[-1]=~s/\.bcf//g;
}if($inputfile=~/\//){
	@parentname=split/\//,$inputfile;
	$parentname[-1]=~s/\.vcf//g;
}elsif($inputfile!~/\//){
	@parentname=split/\//,$abs_inputfilepath;
	$parentname[-1]=~s/\.vcf//g;
}

open IN, "$inputfile";
open OU, ">$parentname[-1]-$progenyname[-1].vcf.loci";
open OUT, ">$parentname[-1]-$progenyname[-1].vcf.loci-seq";
while(<IN>){
        chomp;
        if($_=~/#/ or $_!~/INDEL/){
                next;
        }elsif($_!~/#/){
                my @line=split/\s+/,$_;
                my @vcf_gt=split/\:/,$line[-1];
                my @vcf_allele=split/\//,$vcf_gt[0];
                if($vcf_allele[0]!=$vcf_allele[1]){
                        print OU "$line[0]\t$line[1]\n";
                        print OUT "$line[0]\t$line[1]\t$line[3]\n";
                }else{
                        next;
                }
        }else{
                next;
        }
}
close IN;
close OU;
close OUT;
if(!-e "$progeny_inputfile.csi"){
	`$bcftools index $progeny_inputfile`;
	`$bcftools call -m -f gq -T $parentname[-1]-$progenyname[-1].vcf.loci $progeny_inputfile | awk '\$1!~/^#/' > $progenyname[-1]-$parentname[-1].vcf`;
}if(-e "$progeny_inputfile.csi"){
	`$bcftools call -m -f gq -T $parentname[-1]-$progenyname[-1].vcf.loci $progeny_inputfile | awk '\$1!~/^#/' > $progenyname[-1]-$parentname[-1].vcf`;
}
unlink "$parentname[-1]-$progenyname[-1].vcf.loci";
open IN, "$parentname[-1]-$progenyname[-1].vcf.loci-seq";
while(<IN>){
        chomp;
        my @line=split/\s+/,$_;
        my $l="$line[0]\t$line[1]";
        $hash1{$l}="$_";
}
open IN, "$progenyname[-1]-$parentname[-1].vcf";
while(<IN>){
	chomp;
	if($_=~/INDEL/){
	my @line=split/\s+/,$_;
	my $l="$line[0]\_$line[1]";
	$hash{$l}=1;
	}
}
open IN,"$progenyname[-1]-$parentname[-1].vcf";
open OU,">$progeny_inputfile.vcf.gt";
while(<IN>){
	chomp;
	my @line=split/\s+/,$_;
        my @dp=split /;/,$line[7];
        my @cd=split /:/,$line[-1];
        my @dp1=split /=/,$dp[3];
        my @dp2=split /=/,$dp[4];
        my @dp3=split /,/,$dp2[1];
	my $k="$line[0]\_$line[1]";
	if(exists $hash{$k} && $_!~/INDEL/){ 
		next;
	}else{
		my @b=split/\:/,$line[-1];
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
                }elsif($GT=~/0\/0/ && $line[5]>=20){
                        print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$GT\n";
		}elsif($GT=~/0\/0/ && $line[5]<20){
			print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t_/_\n";
                }elsif($GT=~/\.\/\./){
                        print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t_/_\n";
                }elsif($GT=~/1\/1/ && (($dp3[1]<$HOMOZYGOUS) or ($cd[-1]<$GQ)or ($line[5]<20))){
                        print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t_/_\n";
                }elsif($GT=~/0\/1/ && (($dp3[0]<$HETEROZYGOUS) or ($dp3[1]<$HETEROZYGOUS) or ($cd[-1]<$GQ) or ($line[5]<20))){
                        print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t_/_\n";
                }elsif($GT!~/0\/0/){
                        $b[0]=~s:(\d).*(\d):$hash{$1}/$hash{$2}:g;
                        print OU "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$b[0]\n";
                }
	}
}
close IN;
close OU;
open IN, "$progeny_inputfile.vcf.gt";
open OU, ">$progeny_inputfile.vcf.gt.filter";
while(<IN>){
	chomp;
	my @line=split/\s+/,$_;
	if($_=~/#/){
		next;
	}else{
		my $id="$line[0]\t$line[1]";
		if(exists $hash1{$id} && $line[-1]=~/0\/0/){
			print OU "$hash1{$id}\t.\t$line[-1]\n";
		}else{
			print OU "$_\n";
		}
	}
}
undef %hash;
undef %hash1;
unlink "$parentname[-1]-$progenyname[-1].vcf.loci-seq","$progeny_inputfile.vcf.gt";
rename "$progeny_inputfile.vcf.gt.filter","$progeny_inputfile.vcf.gt";
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
		`$bcftools call -m -r $l1:$l2-$length1 $progeny_inputfile >> $progeny_inputfile\-00`;
	}else{
		next;
	}
}
open IN,"$progeny_inputfile\-00";
open OU,">$progeny_inputfile\-01";
while(<IN>){
	chomp;
	my @line=split/\s+/,$_;
	if($_=~/#/){
		next;
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
	my $id="$line[0]\t$line[1]";
	$hash{$id}=$_;
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
	foreach my $keys (keys %hash){
		my @fd=split/\s+/,$keys;
		if($fd[0]eq$id[0] && $id1[0]<=$fd[1] && $id1[1]>=$fd[1]){
			print OU "$_\n";
		}if($fd[0]eq$id[0] && $id1[0]<=$fd[1] && $id1[1]>$fd[1]){
			print TXT "$_\n";
		}else{
			next;
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
open OU, ">$progenyname[-1]\_$parentname[-1]\.gt";
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
unlink "$progeny_inputfile.homo", "$progeny_inputfile.vcf.gt";

