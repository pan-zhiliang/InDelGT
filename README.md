# InDelGT: an integrated pipeline for extracting InDel genotypes in a hybrid population with next generation DNA sequencing data
# Introduction
InDelGT is used to extract InDel genotypes across a hybrid population with next generation DNA sequencing data for genetic linkage mapping. It takes three steps as following to finish the whole work.  
           1. generating two parental InDel catalogs;  
           2. calling InDel genotypes across a hybrid population;  
           3. analysis of InDel segregation；

# Usage
To run InDelGT, users have to download and install the necessary software packages, including [BWA](https://github.com/lh3/bwa/releases/tag/v0.7.17) and [SAMtools( with BCFtools)](http://samtools.sourceforge.net/), and prepare a parameter setting file, namely `parameters.ini`. The parameter file consists of three parts: folders, parameters and data files. The first part ‘folders’ gives the software paths of InDelGT itself, BWA, SAMtools, BCFtools, the path of Perl program required by InDelGT, the path of fasta file of the reference genome, and the path where the fastq data files of the two parents and all their progeny are saved. The second part ‘parameters’ includes the minimum score of InDel genotypes, the minimum score of the mapping quality of reads, the percentage of missing genotypes required, the reads depth of heterozygote and homozygote，the minimum *P*-value allowed for testing the segregation ratio at a InDel site and the number of threads used for parallel computing. The third part ‘data files’ presents the fasta file of the reference genome and the first and second read files of the two parents and all their progeny. A typical parameter file looks as following:  

    [folders]
    BWA_FOLD:/mnt/sde2/wuhainan/bwa-0.7.17
    SAMTOOLS_FOLD:/mnt/sde2/wuhainan/panzl/a/samtools-1.9
    BCFTOOLS_FOLD:/mnt/sde2/wuhainan/panzl/a/bcftools-1.9
    INDELGT_FOLD:/mnt/sde2/wuhainan/panzl/a/1/InDelGT
    DATAFOLD:/mnt/sde2/wuhainan/panzl/a/1/InDelGT/sample
    [parameters]
    # the minimum score of InDel genotypes
    GQ:30
    # the minimum score of the mapping quality of reads
    MAPQ:5
    # the number of progeny
    PROGENY_NUMBER:20
    # percent of missing genotypes at each InDel
    MISPCT:20
    HETEROZYGOUS_DEPTH:3
    HOMOZYGOUS_DEPTH:5
    PVALUE:0.01
    THREADS:20
    [data files]
    REFERENCE_FILE:reference.fasta
    #When the population is BC population, PARENT1 is the heterozygous parent and PARENT2 is the homozygous parent.
    PARENT1:female.R1.fq  female.R2.fq
    PARENT2:male.R1.fq  male.R2.fq
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

        Usage: perl InDelGT.pl [Options] <type> -o <directory>
        
        Options:
                -p  <type>      the type of population: CP;BC1;BC2;F2 (default: CP)
                -o <directory>  create a directory for storing output file of the InDel genotyping results  
                --help|h        help  
# Test Data

For users to quickly grasp the use of InDelGT, we provide an online test data. All reads data files can be downloaded in a compressed file as [exampledata.tar.gz](https://figshare.com/articles/dataset/sample_tar_gz/15131649), including a reference genome sequence, the first and second reading files of 2 parents and 20 progeny samples. With the parameter file [parameters.ini](https://github.com/pan-zhiliang/InDelGT/blob/main/InDelGT/parameters.ini) attached at this site, we can perform the analysis process by inputting the command: `perl PathToInDelGT/InDelGT.pl`. When the computing finishes, InDels was divided into 7 segregation types: ab×aa, aa×ab, ab×ab, ab×cc, aa×bc, ab×ac or ab×cd. These InDel genotype data for each segregation type are saved in a separate text file.

