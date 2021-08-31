# InDelGT: an integrated pipeline for extracting InDel genotypes in a hybrid population with next generation DNA sequencing data
# Introduction
InDelGT is used to extract InDel genotypes across a hybrid population with next generation DNA sequencing data for genetic linkage mapping. It takes three steps as following to finish the whole work.  
           1. generating two parental InDel catalogs;  
           2. calling InDel genotypes across a hybrid population;  
           3. analysis of InDel segregation；

# Usage
To run InDelGT, users have to download and install the necessary software packages, including [BWA](https://github.com/lh3/bwa/releases/tag/v0.7.17) and [SAMtools( with BCFtools)](http://samtools.sourceforge.net/), and prepare a parameter setting file, namely `parameters.ini`. The parameter file consists of four parts: folders, data files, parameters and fastq files. The first part ‘folders’ gives the software paths of InDelGT itself, BWA, SAMtools, BCFtools, the path of Perl program required by InDelGT, and the path where all the progeny fastq data files are saved. The second part ‘data files’ gives the path of data files such as fasta file of the reference genome, and pair-end sequencing files of the two parents. The third part ‘parameters’ includes the minimum score of InDel genotypes, the minimum score of the mapping quality of reads, the percentage of missing genotypes required, the reads depth of heterozygote and homozygote, the minimum *P*-value allowed for testing the segregation ratio at a InDel site and the number of threads used for parallel computing. The fourth part ‘fastq files’ presents the first and second read files of all progeny. A typical parameter file looks as following:  

    [folders]
    BWA_FOLD:/mnt/sde2/wuhainan/bwa-mem2/bwa-mem2
    SAMTOOLS_FOLD:/mnt/sde2/wuhainan/panzl/a/samtools-1.9/samtools
    BCFTOOLS_FOLD:/mnt/sde2/wuhainan/panzl/a/bcftools-1.9/bcftools
    INDELGT_FOLD:/mnt/sde2/wuhainan/panzl/a/7/CI/InDelGT
    PLSFOLD:/mnt/sde2/wuhainan/panzl/a/7/CI/InDelGT/pls
    PROGENY_FOLD:/home/tong/HuadaData20170417
    [data files]
    REFERENCE_GENOME_FILE:/mnt/sde2/wuhainan/panzl/a/7/CI/Ptrichocarpa_533_v4.0.fa
    MALEPARENT_1_FASTQ_FILE:/mnt/sde2/wuhainan/panzl/a/7/male_merge_1.fq
    MALEPARENT_2_FASTQ_FILE:/mnt/sde2/wuhainan/panzl/a/7/male_merge_2.fq
    FEMALEPARENT_1_FASTQ_FILE:/mnt/sde2/wuhainan/panzl/a/7/female_merge_1.fq
    FEMALEPARENT_2_FASTQ_FILE:/mnt/sde2/wuhainan/panzl/a/7/female_merge_2.fq
    [parameters]
    GQ:30
    MAPQ:10
    PROGENY_NUMBER:117
    MISPCT:10
    HETEROZYGOUS_DEPTH:3
    HOMOZYGOUS_DEPTH:5
    PVALUE:0.01
    THREADS:20
    [fastq files]
    PROGENY1: sample01.R1.fq  sample01.R2.fq
    PROGENY2: sample02.R1.fq  sample02.R2.fq
    PROGENY3: sample03.R1.fq  sample03.R2.fq
    ......
    PROGENY19: sample19.R1.fq  sample19.R2.fq
    PROGENY20: sample20.R1.fq  sample20.R2.fq

  Additionally, users need to install several perl modules, including `Parallel::ForkManager`, `Getopt::Long` and `Statistics::Distributions`. The modules can be easily installed with the module of [cpan](https://www.cpan.org/) if it is installed in advance.  
  When the required software packages are installed and the parameter file is parepared and saved in a work directory, you can go to the work directory and get started with the command:  
  `perl PathToInDelGT/InDelGT.pl -o directory`
    
  Since it will take quite a few hours or even several days to finish a pratical computing, we usually run the command in background as  
  `nohup perl PathToInDelGT/InDelGT.pl -o directory &`

  You can run with the 'help' option (`perl PathToInDelGT/InDelGT.pl -h`) to show the usage of gmRAD:

        Usage: perl InDelGT.pl [Options] -o directory
        
        Options:
                -p  <type>      the type of population: CP;BC1;BC2;F2 (default: CP)
                -o <directory>  create a directory for storing output file of the InDel genotyping results  
                --help|h        help  
# Test Data

For users to quickly grasp the use of InDelGT, we provide an online test data. All reads data files can be downloaded in a compressed file as [exampledata.tar.gz](https://figshare.com/articles/dataset/sample_tar_gz/15131649), including a reference genome sequence, the first and second reading files of 2 parents and 20 progeny samples. With the parameter file [parameters.ini](https://github.com/tongchf/gmRAD/blob/master/parameters.ini) attached at this site, we can perform the analysis process by inputting the command: `perl PathToInDelGT/InDelGT.pl`. When the computing finishes, InDels was divided into 7 segregation types: ab×aa, aa×ab, ab×ab, ab×cc, aa×bc, ab×ac or ab×cd. These InDel genotype data for each segregation type are saved in a separate text file.

