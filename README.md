ADEPT, a dynamic next generation sequencing data error-detection program with trimming 

=======
ADEPT is a program  that dynamically assesses errors within reads by comparing position-specific and neighboring base quality scores with the distribution for the dataset being analyzed. 

-------------
PREREQUISITES
-------------

1. The main program is developed in Perl v 5.8.8.
2. Parallel::ForkManager module from CPAN   
   (http://search.cpan.org/~dlux/Parallel-ForkManager-0.7.9/lib/Parallel/ForkManager.pm)
3. String::Approx module from CPAN   
   (http://search.cpan.org/~jhi/String-Approx-3.27/Approx.pm)
4. R for ploting                 
   (http://www.r-project.org/)                             

-----------
BASIC USAGE
-----------

* Trimming by quality 5 and filtering reads with any ambiguous base or low complexity.

  $   perl ADEPT.pl  -p 'reads1.fastq reads2.fastq' -d out_directory


-----------
Full USAGE
-----------

     perl ADEPT.pl [options] [-u unpaired.fastq] -p 'reads1.fastq reads2.fastq' -d out_directory
     
     Input File: (can use more than once fastq file)
            -u             Unpaired reads
            
            -p             Paired reads in two files and separate by space in quote
     Trim:
            -qE            5" and 3" ends triming # as quality level (0-40) (default 5) for trimming
            -qC            threhold to call a base to be correct (0-1.0) (default = 0.25, higher quality
                           than 25% the nucleotides at that position within the sampled run )
            -qW            threhold to identifying a nucleotide as an error if it falls below a defined 
                           percentage of the quality scores for that position (0-1.0) (default = 0)
            -qMN           ratio of the of the base quality to the qualities of upstream and downstream positions (0-1.0)
                           By default, all qIN ratios must be at least 0.4 to be considered as a potential erroneous base 
                           (i.e. all adjacent qualities must be at least 2.5 times higher than the quality of the position being investigated).
            -qNS           threhold to identify a nucleotide as an potential error if its neighbors' quality falls below a defined 
                           percentage of the quality scores for that neighbors' position within the sampled run (0-1.0) (default = 0.3)

     Filters:
            -min_L         Trimmed sequence length will have at least minimum length (default:50)
            

            					
     Q_Format:
            -ascii         Encoding type: 33 or 64 or autoCheck (default)
                           Type of ASCII encoding: 33 (standard) or 64 (illumina 1.3+)

            -out_ascii     Output encoding. (default: 33)
     Output:
            -prefix        Output file prefix. (default: QC)

            -stats         Statistical numbers output file (default: prefix.stats.txt)

            -d             Output directory.
     Options:
            -t             # of CPUs to run the script (default:2 )

            -split_size    Split the input file into several sub files by sequence number (default: 1000000) 

            -out_non_trim_reads      <bool> Output not trimmed reads to prefix.discard.fastq (default: 0, not output)

            -debug         keep intermediate files


---------------
VERSION HISTORY
---------------

======== Version 1.1
Stable function release.
Features:
- assesses errors within reads by comparing position-specific and neighboring ba  se quality scores with the distribution for the dataset being analyzed. 
- autocheck quality encoding and quality encoding coversion
- multi-threads  (required Parallel::ForkManager)
- input paired end reads aware

-------------
CITATION
-------------
