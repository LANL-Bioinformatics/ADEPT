#!/usr/bin/perl 
######   Required  ########################################################
#          1. Parallel::ForkManager module from CPAN                      #
#          2. String::Approx  module from CPAN                            #
#          3. R for ploting                                               #
###########################################################################
######   Inputs   #########################################################
#          1. Fastq files either paired-end or unpaired reads or both     #
#              Can input multiple library fastq files but only output     #
#              concatenate trimmed fastq files                            #
#          2. Output directory                                            #
#          3. Other options                                               #
###########################################################################
######   Output   #########################################################
#          1.  Two Paired-ends files if input paired-end reads            #
#          2.  One unpaired reads file                                    #
#          3.  statistical text file                                      #
#          4.  quality report pdf file                                    #
###########################################################################
use strict;
use File::Basename;
use Getopt::Long;
use Parallel::ForkManager;
use String::Approx;

my $version=1.1;

my $debug=0;

sub Usage {
print <<"END";
     Usage: perl $0 [options] [-u unpaired.fastq] -p 'reads1.fastq reads2.fastq' -d out_directory
     Version $version
     Input File: (can use more than once)
            -u            Unpaired reads
            
            -p            Paired reads in two files and separate by space in quote
     Trim:
            -qE            5" and 3" ends triming # as quality level (default 5) for trimming
            -qC            threhold to call a base to be correct (default = 0.25, higher quality
                           than 25% the nucleotides at that position within the sampled run )
            -qW            threhold to identifying a nucleotide as an error if it falls below a defined 
                           percentage of the quality scores for that position (default = 0)
            -qMN           ratio of the of the base quality to the qualities of upstream and downstream positions
                           By default, all qIN ratios must be at least 0.4 to be considered as a potential erroneous base 
                           (i.e. all adjacent qualities must be at least 2.5 times higher than the quality of the position being investigated).
            -qNS           threhold to identify a nucleotide as an potential error if its neighbors' quality falls below a defined 
                           percentage of the quality scores for that neighbors' position within the sampled run (default = 0.3)

     Filters:
            -min_L        Trimmed sequence length will have at least minimum length (default:50)
            

            					
     Q_Format:
            -ascii        Encoding type: 33 or 64 or autoCheck (default)
                          Type of ASCII encoding: 33 (standard) or 64 (illumina 1.3+)

            -out_ascii    Output encoding. (default: 33)
     Output:
            -prefix       Output file prefix. (default: QC)

            -stats        Statistical numbers output file (default: prefix.stats.txt)

            -d            Output directory.
     Options:
            -t            # of CPUs to run the script (default:2 )

            -split_size   Split the input file into several sub files by sequence number (default: 1000000) 

            -out_non_trim_reads      <bool> Output not trimmed reads to prefix.discard.fastq (default: 0, not output)

            -debug        keep intermediate files
END
exit;
}

# magic number of quality score
my $highest_illumina_score=41;
my $lowest_illumina_score=0;
# Options Variable initialization
my $thread=2;
my $opt_q=5;
my $opt_min_L=50;
my $minilength=$opt_min_L;
my $opt_avg_cutoff=0;
my $ascii;
my $out_offset=33;
my $subfile_size=1000000;
my $kmer=31; 
my $prefix="QC";
my $plots_file;
my $stats_output;
my $trimmed_reads1_fastq_file;
my $trimmed_reads2_fastq_file;
my $trimmed_unpaired_fastq_file;
my $notrimmed_reads1_fastq_file;
my $notrimmed_reads2_fastq_file;
my @paired_files;
my @unpaired_files;
my $outDir;
my $output_discard;
my $qc_only;
my $nontrimed_output;

#ADEPT OPTIOPN
my $cutlimit=0.3;
my $lowperc=0;
my $upperc=0.25;
 my $discperc=0.35;
my $cutlimit1= $cutlimit;
my $cutlimit2= $cutlimit;


# Options
GetOptions("qE=i"         => \$opt_q,
           "qC=i"         => \$upperc,
           "qW=i"         => \$lowperc,
           "qMN=i"        => \$discperc, 
           "qNS=i"        => \$cutlimit, 
           "min_L=i"      => \$opt_min_L,
           "avg_q=f"      => \$opt_avg_cutoff,
           "p=s"          => \@paired_files,
           "u=s"          => \@unpaired_files,
           "ascii=i"      => \$ascii,
           "out_ascii=i"  => \$out_offset,
           "t|threads=i"  => \$thread,
           'split_size=i' => \$subfile_size,
           'prefix=s'     => \$prefix,
           'd=s'          => \$outDir,
           'stats=s'      => \$stats_output,
           'out_non_trim_reads'      => \$nontrimed_output,
           'debug'        => \$debug,
           "help|?"       => sub{Usage()} );

die "Missing input files.\n", &Usage() unless @unpaired_files or @paired_files;
die "Missing output directory.\n",&Usage() unless $outDir;

###### Output file initialization #####
# temp files for plotting
my $quality_matrix="$outDir/$prefix.quality.matrix";
my $avg_quality_histogram="$outDir/$prefix.for_qual_histogram.txt";
my $base_matrix="$outDir/$prefix.base.matrix";
my $nuc_composition_file="$outDir/$prefix.base_content.txt";
my $length_histogram="$outDir/$prefix.length_count.txt";

# output files
$plots_file="$outDir/${prefix}_qc_report.pdf" if (!$plots_file);
$trimmed_unpaired_fastq_file="$outDir/$prefix.unpaired.trimmed.fastq" if (!$trimmed_unpaired_fastq_file);
$trimmed_reads1_fastq_file="$outDir/$prefix.1.trimmed.fastq" if (!$trimmed_reads1_fastq_file);
$trimmed_reads2_fastq_file="$outDir/$prefix.2.trimmed.fastq" if (!$trimmed_reads2_fastq_file);
$notrimmed_reads1_fastq_file="$outDir/$prefix.1.qc.fastq" if (!$trimmed_reads1_fastq_file);
$notrimmed_reads2_fastq_file="$outDir/$prefix.2.qc.fastq" if (!$trimmed_reads2_fastq_file);


$stats_output="$outDir/$prefix.stats.txt" if (!$stats_output);
   
#######################################

######   Output check  ################
if (! -e $outDir)
{
  mkdir $outDir;
}

if (-e $trimmed_reads1_fastq_file)
{
   print "The output $trimmed_reads1_fastq_file file exists and will be overwritten.\n";
   system ("rm $trimmed_reads1_fastq_file");
}
if (-e $trimmed_reads2_fastq_file)
{
   print "The output $trimmed_reads2_fastq_file file exists and will be overwritten.\n";
   system ("rm $trimmed_reads2_fastq_file");
}
if ($notrimmed_reads1_fastq_file)
{
if (-e $notrimmed_reads1_fastq_file)
{
   print "The output $notrimmed_reads1_fastq_file file exists and will be overwritten.\n";
   system ("rm $notrimmed_reads1_fastq_file");
}
}
if ( $notrimmed_reads2_fastq_file)
{
if (-e $notrimmed_reads2_fastq_file)
{
   print "The output $notrimmed_reads2_fastq_file file exists and will be overwritten.\n";
   system ("rm $notrimmed_reads2_fastq_file");
}
}
if (-e $trimmed_unpaired_fastq_file)
{
   print "The output $trimmed_unpaired_fastq_file file exists and will be overwritten.\n";
   system ("rm $trimmed_unpaired_fastq_file");
}
if (-e $plots_file)
{
   print "The output $plots_file file exists and will be overwritten.\n";
   system ("rm $plots_file");
}
#######################################



my ( $total_count,$total_count_1, $total_count_2, $total_num, $trimmed_num,$total_raw_seq_len, $total_trimmed_seq_len);
my ( $trim_seq_len_std, $trim_seq_len_avg, $max, $min, $median, $numOfReadsWithAdapter);
my ( $paired_seq_num, $total_paired_bases );
my ( @split_files, @split_files_2);
my (%error);
my %position;
my %AverageQ;
my %base_position;
my %base_content;
my %len_hash;
my ( $i_file_name, $i_path, $i_suffix );
  
  foreach my $input (@unpaired_files,@paired_files){
     print "Processing $input file\n";
     #print $STATS_fh "Processing $input file\n";
     my ($reads1_file,$reads2_file) = split /\s+/,$input;
  
     # check file 
     if(&file_check($reads1_file)<0) { die "The file $reads1_file doesn't exist or empty.\n";}
     if(&file_check($reads2_file)<0 and $reads2_file) { die "The file $reads2_file doesn't exist or empty.\n";}

     # check quality offset
     if (! $ascii){$ascii = &checkQualityFormat($reads1_file)}
    #baseline quality

    my $basequality=&build_initial_quality_matrix($reads1_file);
    #split
    ($total_count_1,@split_files) = &split_fastq($reads1_file,$outDir,$subfile_size);
    ($total_count_2,@split_files_2) = &split_fastq($reads2_file,$outDir,$subfile_size) if ($reads2_file);
     $total_count += $total_count_1 + $total_count_2;

    my $pm = new Parallel::ForkManager($thread);

    $pm -> run_on_finish ( # called BEFORE the first call to start()
      sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $nums_ref) = @_;
      
        # retrieve data structure from child
        if (defined($nums_ref)) {  # children are not forced to send anything
          my ($total_score, $avg_score);
          $total_num += $nums_ref->{raw_seq_num};
          $trimmed_num += $nums_ref->{trim_seq_num};
          $total_raw_seq_len += $nums_ref->{total_raw_seq_len};
          $total_trimmed_seq_len += $nums_ref->{total_trim_seq_len};
          my $processed_num= $nums_ref->{raw_seq_num};
          $trim_seq_len_std = $nums_ref->{trim_seq_len_std};
          $trim_seq_len_avg = $nums_ref->{trim_seq_len_avg};
          $max = $nums_ref->{max};
          $min = $nums_ref->{min};
          $median = $nums_ref->{median};
          $paired_seq_num += $nums_ref->{paired_seq_num};
          $total_paired_bases +=  $nums_ref->{total_paired_bases};
          $numOfReadsWithAdapter += $nums_ref->{numOfReadsWithAdapter};
          my %temp_avgQ = %{$nums_ref->{ReadAvgQ}};
          map {$AverageQ{$_}->{bases} += $temp_avgQ{$_}->{basesNum};
               $AverageQ{$_}->{reads} += $temp_avgQ{$_}->{readsNum}; 
              } keys %temp_avgQ; 
          my %temp_position= %{$nums_ref->{qual}};
          my %temp_base_position= %{$nums_ref->{base}};
          foreach my $pos (1..($total_raw_seq_len/$total_num))
          {
                for my $score ($lowest_illumina_score..$highest_illumina_score)
                { 
                    $position{$pos}->{$score} += $temp_position{$pos}->{$score};
                    $total_score +=  $score * $temp_position{$pos}->{$score}; 
                }
                for my $nuc ("A","T","C","G","N")
                {
                    $base_position{$pos}->{$nuc} += $temp_base_position{$pos}->{$nuc};
                }
          }
          my %tmp_base_content = %{$nums_ref->{Base_content}};
          for my $nuc ("A","T","C","G","N","GC")
          {
              while (my ($key, $value)= each %{$tmp_base_content{$nuc}})
              {
                 $base_content{$nuc}->{$key} += $value;
              }
          } 
          while (my ($key, $value)= each %{$nums_ref->{ReadLen}} )
          {
		         $len_hash{$key} += $value;   
          }

          #print $STATS_fh " Processed $total_num/$total_count\n";
          print "Processed $total_num/$total_count\n";
          printf (" Post Trimming Length(Mean, Std, Median, Max, Min) of %d reads with Overall quality %.2f\n",$processed_num, $total_score/$nums_ref->{total_trim_seq_len});
          printf (" (%.2f, %.2f, %.1f, %d, %d)\n",$trim_seq_len_avg,$trim_seq_len_std,$median,$max,$min);
          #unlink $split_files[$ident];
          #unlink $split_files_2[$ident] if ($split_files_2[$ident]);
         
        } else {  # problems occuring during storage or retrieval will throw a warning
          print qq|No message received from child process $pid! on $ident\n|;
        }
      }
    );

  foreach my $i(0..$#split_files)
  {
    $pm->start($i) and next;
    my $hash_ref = &qc_process($split_files[$i],$split_files_2[$i],$basequality);
    $pm->finish(0, $hash_ref);
   
  }
  $pm->wait_all_children;
  # clean up
  foreach my $i(0..$#split_files)
  {
    unlink $split_files[$i];
    unlink $split_files_2[$i] if ($split_files_2[$i]);
  }
} #end foreach $input


# concatenate each thread's trimmed reads and files clean up.
if (! $qc_only)
{
    if (@unpaired_files) 
    {
      foreach my $input (@unpaired_files){
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$input", qr/\.[^.]*/ );
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$i_file_name", qr/\.[^.]*/ ) if ($i_suffix =~ /gz$/);
    
        if (system("cat $outDir/${i_file_name}_?????_trim.fastq >> $trimmed_unpaired_fastq_file")) { die "cat failed: $!" }
        `rm $outDir/${i_file_name}_?????_trim.fastq`;
      }
    }
    
    if (@paired_files)
    {
      foreach my $input (@paired_files){
        my ($reads1_file,$reads2_file) = split /\s+/,$input;
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$reads1_file", qr/\.[^.]*/ );
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$i_file_name", qr/\.[^.]*/ ) if ($i_suffix =~ /gz$/);
        if (system("cat $outDir/${i_file_name}_?????_trim.fastq >> $trimmed_reads1_fastq_file")) { die "cat failed: $!" }
        `rm $outDir/${i_file_name}_?????_trim.fastq`;
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$reads2_file", qr/\.[^.]*/ );
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$i_file_name", qr/\.[^.]*/ ) if ($i_suffix =~ /gz$/);
        if (system("cat $outDir/${i_file_name}_?????_trim.fastq >> $trimmed_reads2_fastq_file")) { die "cat failed: $!" }
        `rm $outDir/${i_file_name}_?????_trim.fastq`;
        if (system("cat $outDir/${i_file_name}_?????_trim_unpaired.fastq >> $trimmed_unpaired_fastq_file")) { die "cat failed: $!" }
        `rm $outDir/${i_file_name}_?????_trim_unpaired.fastq`;
      }
    } 

  if ($nontrimed_output) 
  {
 

    if (@paired_files)
    {
      foreach my $input (@paired_files){
        my ($reads1_file,$reads2_file) = split /\s+/,$input;
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$reads1_file", qr/\.[^.]*/ );
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$i_file_name", qr/\.[^.]*/ ) if ($i_suffix =~ /gz$/);
        if (system("cat $outDir/${i_file_name}_?????_qc.fastq >> $notrimmed_reads1_fastq_file")) { die "cat failed: $!" }
        `rm $outDir/${i_file_name}_?????_qc.fastq`;
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$reads2_file", qr/\.[^.]*/ );
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$i_file_name", qr/\.[^.]*/ ) if ($i_suffix =~ /gz$/);
        if (system("cat $outDir/${i_file_name}_?????_qc.fastq >> $trimmed_reads2_fastq_file")) { die "cat failed: $!" }
        `rm $outDir/${i_file_name}_?????_qc.fastq`;
      }
    }

   
  }
}

&print_quality_report_files();
&print_final_stats();
&plot_by_R();
unless ($debug){
unlink $nuc_composition_file;
unlink $quality_matrix;
unlink $base_matrix;
unlink $avg_quality_histogram;
unlink $length_histogram;
}
exit(0);

sub build_initial_quality_matrix 
{
  print "start build_initial_quality_matrix\n";
    my $fastqFile=shift;
    open (my $fh, "$fastqFile") or die "$fastqFile $!";
    my %basequal;
    my $total_reads=0;
    while(<$fh>)
    {
        last if ($total_reads>1000000);
        my $name=$_;
        my $seq=<$fh>;
        $seq =~ s/\n//g;
        while ($seq !~ /\+/)
        {
           $seq .= <$fh>;
           $seq =~ s/\n//g;
        }
        my $q_id_pos=index($seq,"+");
        $seq = substr($seq, 0, $q_id_pos);
        my $seq_len = length $seq;
        my $qual_seq=<$fh>;
        $qual_seq =~ s/\n//g;
        my $qual_seq_len = length $qual_seq;
        while ( $qual_seq_len < $seq_len )
        {
           last if ( $qual_seq_len == $seq_len);
           $qual_seq .= <$fh>;
           $qual_seq =~ s/\n//g;
           $qual_seq_len = length $qual_seq;
        }
        my @qual_seq=split //, $qual_seq;

        if (rand() <=0.2) {
           $total_reads++;
           for my $pos(0..$#qual_seq)
           {
              push @{$basequal{$pos}}, ord($qual_seq[$pos])-$ascii;
           }
        }
    }
    close $fh;
    @{$basequal{$_}}= sort {$a <=> $b} @{$basequal{$_}} foreach (keys %basequal);

      print "done build_initial_quality_matrix\n";
    return \%basequal;
}



sub print_final_stats{
    open (my $fh, ">$stats_output") or die "$!\t$stats_output\n";
    # stats
    print $fh "\nBefore Trimming\n";
    print $fh "Reads #:\t$total_num\n";
    print $fh "Total bases:\t$total_raw_seq_len\n";
    printf $fh ("Reads Length:\t%.2f\n",$total_raw_seq_len/$total_num);
    
    print $fh "\nAfter Trimming\n";
    printf $fh ("Reads #:\t\%d (%.2f %%)\n",$trimmed_num, $trimmed_num/$total_num*100);
    printf $fh ("Total bases:\t\%d (%.2f %%)\n",$total_trimmed_seq_len,$total_trimmed_seq_len/$total_raw_seq_len*100);
    if ($trimmed_num)
    {
      printf $fh ("Mean Reads Length:\t%.2f\n",$total_trimmed_seq_len/$trimmed_num); 
    }
    else
    {
      printf $fh "Mean Reads Length:\t0\n";
    }
    
    if (@paired_files){
      printf $fh ("  Paired Reads #:\t\%d (%.2f %%)\n",$paired_seq_num, $paired_seq_num/$trimmed_num*100);
      printf $fh ("  Paired total bases:\t\%d (%.2f %%)\n",$total_paired_bases,$total_paired_bases/$total_trimmed_seq_len*100);
      printf $fh ("  Unpaired Reads #:\t\%d (%.2f %%)\n", $trimmed_num - $paired_seq_num, ($trimmed_num - $paired_seq_num)/$trimmed_num*100);
      printf $fh ("  Unpaired total bases:\t\%d (%.2f %%)\n", $total_trimmed_seq_len - $total_paired_bases , ($total_trimmed_seq_len - $total_paired_bases)/$total_trimmed_seq_len*100);
    }
    
    printf $fh ("\nDiscarded reads #:\t\%d (%.2f %%)\n", $total_num - $trimmed_num , ($total_num - $trimmed_num)/$total_num*100);
    printf $fh ("Trimmed bases:\t\%d (%.2f %%)\n", $total_raw_seq_len - $total_trimmed_seq_len, ($total_raw_seq_len - $total_trimmed_seq_len)/$total_raw_seq_len*100);
    
    return (0);
}


sub plot_by_R
{
    open (R,">$outDir/tmp$$.R");
    print R <<RSCRIPT; 
def.par <- par(no.readonly = TRUE) # get default parameters

pdf(file = \"$plots_file\",width = 10, height = 8)
#trimmed summary
txt<-c(
"Before Trimming",
paste("Reads #:", prettyNum(\"$total_num\",big.mark=",")),
paste("Total bases:", prettyNum(\"$total_raw_seq_len\",big.mark=",")),
paste("Reads Length:",sprintf("%.2f",$total_raw_seq_len/$total_num)),
" ",
"After Trimming",
paste("Reads #:", prettyNum(\"$trimmed_num\",big.mark=","),sprintf("(%.2f %%)", $trimmed_num/$total_num*100)),
paste("Total bases:",prettyNum(\"$total_trimmed_seq_len\",big.mark=","),sprintf("(%.2f %%)",$total_trimmed_seq_len/$total_raw_seq_len*100))
)



if($trimmed_num>0)
{
txt<-c(txt,
paste("Mean Reads Length:", sprintf("%.2f",$total_trimmed_seq_len/$trimmed_num))
)
}

if($paired_seq_num >0){
     txt <- c(txt, 
 paste("  Paired Reads #: ", prettyNum(\"$paired_seq_num\",big.mark=","), sprintf("(%.2f %%)", $paired_seq_num/$trimmed_num*100)),
 paste("  Paired total bases: ", prettyNum(\"$total_paired_bases\",big.mark=","),sprintf("(%.2f %%)",$total_paired_bases/$total_trimmed_seq_len*100)),
 paste("  Unpaired Reads: ", prettyNum($trimmed_num - $paired_seq_num,big.mark=","),sprintf("(%.2f %%)", ($trimmed_num - $paired_seq_num)/$trimmed_num*100)),
 paste("  Unpaired total bases: ",prettyNum($total_trimmed_seq_len - $total_paired_bases,big.mark=","),sprintf("(%.2f %%)" , ($total_trimmed_seq_len - $total_paired_bases)/$total_trimmed_seq_len*100))
 )
}

txt<-  c( txt, " ")
    txt<-  c( txt,
   paste("Discarded reads #: ", prettyNum($total_num - $trimmed_num,big.mark=","), sprintf("(%.2f %%)", ($total_num - $trimmed_num)/$total_num*100)),
   paste("Discarded bases: ",prettyNum($total_raw_seq_len - $total_trimmed_seq_len,big.mark=","),sprintf("(%.2f %%)",  ($total_raw_seq_len - $total_trimmed_seq_len)/$total_raw_seq_len*100))
  )


plot(1:80,xaxt=\'n\',yaxt=\'n\',type=\'n\',ylab=\'\',xlab=\'\')
#title(paste(\"$prefix\",\"QC report\"),sub = 'DOE Joint Genome Institute/Los Alamos National Laboratory', adj = 0.5, col.sub='darkblue',font.sub=2,cex.sub=0.8)
title("QC stats")
for (i in seq(along=txt))
{
 text(1,(85-i*5),txt[i],adj=0,font=2)
}
 

#lenght histogram
lengthfile<-read.table(file=\"$length_histogram\")
lengthList<-as.numeric(lengthfile\$V1)
lengthCount<-as.numeric(lengthfile\$V2)
lenAvg<-sum(lengthList * lengthCount)/sum(lengthCount)
lenStd<-sqrt(sum(((lengthList - lenAvg)**2)*lengthCount)/sum(lengthCount))
lenMax<-max(lengthList[lengthCount>0])
lenMin<-min(lengthList[lengthCount>0])
barplot(lengthCount/1000000,names.arg=lengthList,xlab=\"Length\",ylab=\"Count (millions)\",main=\"Reads Length Histogram\",cex.names=0.8)
legend.txt<-c(paste(\"Mean\",sprintf (\"%.2f\",lenAvg),\"±\",sprintf (\"%.2f\",lenStd)),paste(\"Max\",lenMax),paste(\"Min\",lenMin))
legend('topleft',legend.txt,bty='n')

#readGC plot
baseP<-read.table(file=\"$nuc_composition_file\")
Apercent<-baseP\$V2[which(baseP\$V1==\"A\")]
ApercentCount<-baseP\$V3[which(baseP\$V1==\"A\")]
Tpercent<-baseP\$V2[which(baseP\$V1==\"T\")]
TpercentCount<-baseP\$V3[which(baseP\$V1==\"T\")]
Cpercent<-baseP\$V2[which(baseP\$V1==\"C\")]
CpercentCount<-baseP\$V3[which(baseP\$V1==\"C\")]
Gpercent<-baseP\$V2[which(baseP\$V1==\"G\")]
GpercentCount<-baseP\$V3[which(baseP\$V1==\"G\")]
#Npercent<-baseP\$V2[which(baseP\$V1==\"N\")]
#NpercentCount<-baseP\$V3[which(baseP\$V1==\"N\")]
GCpercent<-baseP\$V2[which(baseP\$V1==\"GC\")]
GCpercentCount<-baseP\$V3[which(baseP\$V1==\"GC\")]
aAvg<-sum(Apercent * ApercentCount)/sum(ApercentCount)
aStd<-sqrt(sum(((Apercent - aAvg)**2)*ApercentCount)/sum(ApercentCount))
tAvg<-sum(Tpercent * TpercentCount)/sum(TpercentCount)
tStd<-sqrt(sum(((Tpercent - tAvg)**2)*TpercentCount)/sum(TpercentCount))
cAvg<-sum(Cpercent * CpercentCount)/sum(CpercentCount)
cStd<-sqrt(sum(((Cpercent - cAvg)**2)*CpercentCount)/sum(CpercentCount))
gAvg<-sum(Gpercent * GpercentCount)/sum(GpercentCount)
gStd<-sqrt(sum(((Gpercent - gAvg)**2)*GpercentCount)/sum(GpercentCount))
#nAvg<-sum(Npercent * NpercentCount)/sum(NpercentCount)
#nStd<-sqrt(sum(((Npercent - nAvg)**2)*NpercentCount)/sum(NpercentCount))
gcAvg<-sum(GCpercent * GCpercentCount)/sum(GCpercentCount)
gcStd<-sqrt(sum(((GCpercent - gcAvg)**2)*GCpercentCount)/sum(GCpercentCount))
GCaggregate<-tapply(GCpercentCount,list(cut(GCpercent,breaks=c(seq(0,100,1)))),FUN=sum)
Aaggregate<-tapply(ApercentCount,list(cut(Apercent,breaks=c(seq(0,100,1)))),FUN=sum)
Taggregate<-tapply(TpercentCount,list(cut(Tpercent,breaks=c(seq(0,100,1)))),FUN=sum)
Caggregate<-tapply(CpercentCount,list(cut(Cpercent,breaks=c(seq(0,100,1)))),FUN=sum)
Gaggregate<-tapply(GpercentCount,list(cut(Gpercent,breaks=c(seq(0,100,1)))),FUN=sum)


par(fig=c(0,0.75,0,1),mar=c(5,4,4,2),xpd=FALSE,cex.main=1.2)
plot(GCaggregate/1000000,xlim=c(0,100),type=\"h\",lwd=4, main=\"Reads GC content\",xlab=\"GC (%)\",ylab=\"Number of reads (millions)\",lend=2)
legend.txt<-c(paste(\"GC\",sprintf (\"%.2f%%\",gcAvg),\"±\",sprintf (\"%.2f\",gcStd)))
legend('topright',legend.txt,bty='n')

par(fig=c(0.75,1,0.75,1), mar=c(3, 2, 2, 2),new=TRUE,cex.main=1)
legend.txt<-c(paste(\"A\",sprintf (\"%.2f%%\",aAvg),\"±\",sprintf (\"%.2f\",aStd)))
plot(Aaggregate/1000000,xlim=c(0,50),type=\"h\",lwd=2,main=legend.txt,xlab="",ylab="",)

par(fig=c(0.75,1,0.5,0.75),mar=c(3, 2, 2, 2),new=TRUE)
legend.txt<-c(paste(\"T\",sprintf (\"%.2f%%\",tAvg),\"±\",sprintf (\"%.2f\",tStd)))
plot(Taggregate/1000000,xlim=c(0,50),type=\"h\",lwd=2,main=legend.txt,xlab="",ylab="")

par(fig=c(0.75,1,0.25,0.5),mar=c(3, 2, 2, 2),new=TRUE)
legend.txt<-c(paste(\"C\",sprintf (\"%.2f%%\",cAvg),\"±\",sprintf (\"%.2f\",cStd)))
plot(Caggregate/1000000,xlim=c(0,50),type=\"h\",lwd=2,main=legend.txt,xlab="",ylab="")

par(fig=c(0.75,1,0,0.25),mar=c(3, 2, 2, 2),new=TRUE)
legend.txt<-c(paste(\"G\",sprintf (\"%.2f%%\",gAvg),\"±\",sprintf (\"%.2f\",gStd)))
plot(Gaggregate/1000000,xlim=c(0,50),type=\"h\",lwd=2,main=legend.txt,xlab="",ylab="")
par(def.par)#- reset to default


#ATCG composition per base ATCG plot
baseM<-read.table(file=\"$base_matrix\")
aBase<-baseM\$V1
tBase<-baseM\$V2
cBase<-baseM\$V3
gBase<-baseM\$V4
nBase<-baseM\$V5

aPer<-(aBase/rowSums(baseM))*100
tPer<-(tBase/rowSums(baseM))*100
cPer<-(cBase/rowSums(baseM))*100
gPer<-(gBase/rowSums(baseM))*100

xpos<-seq(1,length(aBase),1)
plot(xpos,aPer,col=\'green3\',type=\'l\',xaxt=\'n\',xlab=\'Base\',ylab=\'Base content (%)\',ylim=c(1,100))
lines(xpos,tPer,col=\'red\')
lines(xpos,cPer,col=\'blue\')
lines(xpos,gPer,col=\'black\')
axis(1,at=xpos,labels=xpos)
legend(\'topright\',c(\'A\',\'T\',\'C\',\'G\'),col=c(\'green3\',\'red\',\'blue\',\'black\'),box.col=0,lwd=1)
title(\"Nucleotide Content Per Cycle\")




# read avg quality count barplot 
Qhist_file<-read.table(file=\"$avg_quality_histogram\",header=TRUE)
par(mar=c(5,4,4,4))
cumulate<-cumsum(Qhist_file\$readsNum)
plot(Qhist_file\$Score,Qhist_file\$readsNum/1000000,type=\'h\',xlim=c(max(Qhist_file\$Score),min(Qhist_file\$Score)),xlab=\"Avg Score\", ylab=\"Reads Number (millions)\",lwd=12,lend=2)
title('Reads Average Quality Histogram')
par(new=TRUE)
plot(Qhist_file\$Score,cumulate/sum(Qhist_file\$readsNum)*100,type='l',xlim=c(max(Qhist_file\$Score),min(Qhist_file\$Score)),yaxt='n',xaxt='n',ylab="",xlab="",col='blue',lwd=3)
axis(4,col='blue',col.ticks='blue',col.axis='blue')
mtext(side=4,'Cumulative Percentage',line=2,col='blue')
Qover20Reads<-sum(as.numeric(Qhist_file\$readsNum[Qhist_file\$Score>=20]))
Qover20ReadsPer<-sprintf(\"%.2f%%\",Qover20Reads/sum(Qhist_file\$readsNum)*100)
Qover20Bases<-sum(as.numeric(Qhist_file\$readsBases[Qhist_file\$Score>=20]))
Qover20AvgLen<-sprintf(\"%.2f\",Qover20Bases/Qover20Reads)
mtext(side=3,paste(\"Number of Q>=20 reads:\",prettyNum(Qover20Reads,big.mark=\",\"),\"(\",Qover20ReadsPer,\")\",\", mean Length:\",Qover20AvgLen),adj=0,cex=0.8,line=0.5)
par(def.par)#- reset to default


# read in matrix file for the following three plots
z<-as.matrix(read.table(file=\"$quality_matrix\"));
x<-1:nrow(z)
y<-1:ncol(z)
y<-y-1

#quality boxplot per base
is.wholenumber <- function(x, tol = .Machine\$double.eps^0.5)  abs(x - round(x)) < tol
plot(1:length(x),x,type=\'n\',xlab=\"Position\",ylab=\"Quality score\", ylim=c(0,max(y)+1),xaxt=\'n\')
axis(1,at=x,labels=TRUE)
title(\"Quality Boxplot Per Cycle\")

for (i in 1:length(x)) {
  total<-sum(z[i,])
  qAvg<-sum(y*z[i,])/total
  if (is.wholenumber(total/2))
  {
     med<-( min(y[cumsum((z[i,]))>=total/2]) + min(y[cumsum((z[i,]))>=total/2+1]) )/2
  }
  else
  {
     med<-min(y[cumsum((z[i,]))>=ceiling(total/2)])
  }

  if (is.wholenumber(total/4))
  {
     Q1<-( min(y[cumsum((z[i,]))>=total/4]) + min(y[cumsum((z[i,]))>=total/4+1]) )/2
  }
  else
  {
     Q1<-min(y[cumsum((z[i,]))>=round(total/4)])
  }

  if (is.wholenumber(total/4*3))
  {
     Q3<-( min(y[cumsum((z[i,]))>=total/4*3]) + min(y[cumsum((z[i,]))>=total/4*3+1]) )/2
  }
  else
  {
     Q3<-min(y[cumsum((z[i,]))>=round(total/4*3)])
  }
  maxi<-max(y[z[i,]>0])
  mini<-min(y[z[i,]>0])
  #if (Q1 == 'Inf') {Q1 = maxi}
  if (Q3 == \'Inf\') {Q3 = maxi}
  IntQ<-Q3-Q1
  mini<-max(mini,Q1-1.5*IntQ)
  maxi<-min(maxi,Q3+1.5*IntQ)
  rect(i-0.4,Q1,i+0.4,Q3,col=\'bisque\')
  lines(c(i,i),c(Q3,maxi),lty=2)
  lines(c(i,i),c(mini,Q1),lty=2)
  lines(c(i-0.4,i+0.4),c(mini,mini))
  lines(c(i-0.4,i+0.4),c(maxi,maxi))
  lines(c(i-0.4,i+0.4),c(med,med))
  #points(i,qAvg,col=\'red\')
  reads_num<-prettyNum($trimmed_num,big.mark=",")
  reads_base<-prettyNum($total_trimmed_seq_len,big.mark=",")
  abline(h=20, col = \"gray60\")
  legend(\"bottomleft\",c(paste(\"# Reads: \",reads_num),paste(\"# Bases:\",reads_base)))
## for outliers
#points()
}

#quality 3D plot
persp(x,y,z/1000000,theta = 50, phi = 30, expand = 0.7, col = \"#0000ff22\",ntick=10,ticktype=\"detailed\",xlab=\'Position\',ylab=\'Score\',zlab=\"\",r=6,shade=0.75)
mtext(side=2, \"Frequency (millions)\",line=2)
title(\"Quality 3D plot. (Position vs. Score vs. Frequency)\")

#Quality count bar plot
col<-colSums(z)
less30columnNum<-length(col)-$highest_illumina_score+30-1
atleast30columnNum<-$highest_illumina_score-30+1
color<-c(rep(\'blue\',less30columnNum),rep(\'darkgreen\',atleast30columnNum))
over30per<-sprintf(\"%.2f%%\",sum(col[(less30columnNum+1):length(col)])/sum(col)*100)
countInM<-col/1000000
plot(seq($lowest_illumina_score,$highest_illumina_score,1),countInM,col=color,type='h',ylab=\"Total (million)\",xlab=\"Q score\",lwd=12,lend=2,bty='n')
abline(v=29.5,col='darkgreen')
text(30,(max(countInM)-min(countInM))*0.9,labels=\">=Q30\",cex=0.8,adj=0,col=\'darkgreen\')
text(30,(max(countInM)-min(countInM))*0.85,labels=over30per,cex=0.8,adj=0,col=\'darkgreen\')
title(\"Quality report\")

tmp<-dev.off()

quit()

RSCRIPT

    close R;
    system ("R --vanilla --silent --quiet < $outDir/tmp$$.R 1>/dev/null");
    unless ($debug){system ("rm $outDir/tmp$$.*");}
    system ("rm $outDir/Rplots.pdf") if (-e "$outDir/Rplots.pdf");
    system ("rm Rplots.pdf") if (-e "Rplots.pdf");
}

sub print_quality_report_files{
    open (OUT, ">$quality_matrix");
    open (BASE, ">$base_matrix");
    foreach my $pos2(sort {$a<=>$b} keys %position){
        my $q_string;
        my $n_string;
        for ($lowest_illumina_score..$highest_illumina_score)
        {
            if ($position{$pos2}->{$_}){
              $q_string .=  $position{$pos2}->{$_}."\t";
            }
            else
            {
              $q_string .= "0\t";
            }
        }
        for my $base ("A","T","C","G","N")
        { 
            if ($base_position{$pos2}->{$base}){
              $n_string .= $base_position{$pos2}->{$base}."\t";
            }
            else
            {
              $n_string .= "0\t";
            }
        }
        $q_string =~ s/\t$/\n/;
        $n_string =~ s/\t$/\n/;
        print OUT $q_string;
        print BASE $n_string;
    }
    close OUT;
    close BASE;

    open (OUT2, ">$avg_quality_histogram");
    my $Key=$highest_illumina_score;
    my $h_print_string;
    while ($Key >= 0)
    {
        if(defined $AverageQ{$Key}->{reads}){
           $h_print_string .= "$Key\t".
           $AverageQ{$Key}->{reads}. "\t".$AverageQ{$Key}->{bases}."\n";
        }else{
           $h_print_string .= "$Key\t".
              "0\t".
              "0\n";
        }
        --$Key;
    }
    print OUT2 "Score\treadsNum\treadsBases\n";
    print OUT2 $h_print_string; 
    close OUT2;

    open (OUT3, ">$nuc_composition_file");
    for my $nuc ("A","T","C","G","N","GC")
    {
      foreach my $key ( sort {$a<=>$b} keys %{$base_content{$nuc}})
      {
             print OUT3 "$nuc\t$key\t${$base_content{$nuc}}{$key}\n";
      }
    } 
    close OUT3;
 
    open (OUT4,">$length_histogram");
	my @len_list= sort {$a<=>$b} keys %len_hash;
    for my $key (1..$len_list[-1])
    {
	    if ($len_hash{$key}) 
		{
            print OUT4 "$key\t$len_hash{$key}\n";
		}
		else
		{
			print OUT4 "$key\t0\n";
		}
    }
    close OUT4;
}

sub file_check 
{
    #check file exist and non zero size
    my $file=shift;
    my $exist=-1;
    if (-e $file) {$exist=1};
    if (-z $file) {$exist=-1};
    return $exist;
}

sub qc_process {
  my ($input, $input2, $basequality) = @_;
  my ($h1,$s,$s_trimmed,$h2,$q, $q_trimmed); my $len=0; my $trim_len=0;
  my ($r2_h1,$r2_s,$r2_s_trimmed,$r2_h2,$r2_q,$r2_q_trimmed); my $r2_len=0; my $r2_trim_len=0;
  my ($q1_nontrim, $s1_nontrim,$q2_nontrim,$s2_nontrim);
  my %stats;
  my $seq_r;
  my $avg_q;
  my ($pos5,$pos3);
  my $numOfReadsWithAdapter;
  my ($raw_seq_num,$total_raw_seq_len);
  my ($trim_seq_num_1,$trim_seq_num_2,$trim_seq_len, $total_trim_seq_len,@trim_seq_len);
  my ($paired_seq_num,$total_paired_bases); 
  my (%tmp1,%tmp2);
  my ($drop_1,$drop_2)=(1,1);
  my ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$input", qr/\.[^.]*/ );
  my $trim_output_1="$outDir/${i_file_name}_trim.fastq";
   my $qc_output_1="$outDir/${i_file_name}_qc.fastq";

  my $qc_output_2;
  my $trim_output_2;
  my $trim_output_unpaired;
  my $trim_output_discard="$outDir/${i_file_name}_trim_discard.fastq";
  open(IN,"$input") or die "$input\t$!";
 if ($nontrimed_output) { open(OUTQC,"> $qc_output_1"); }
 if (! $qc_only)
  {
     open(OUT,"> $trim_output_1");
     open(DISCARD,">$trim_output_discard") if ($output_discard);
  }
  if ($input2) # paired mate
  {
     open(IN2,"$input2") or die "$input2\t$!";
     my ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$input2", qr/\.[^.]*/ );
     $trim_output_2="$outDir/${i_file_name}_trim.fastq";
     my $qc_output_2="$outDir/${i_file_name}_qc.fastq";
     $trim_output_unpaired="$outDir/${i_file_name}_trim_unpaired.fastq";
   
  if ($nontrimed_output)  { open(OUTQC2,"> $qc_output_2"); }

  
   if (! $qc_only)
     {
         open(OUT2,"> $trim_output_2");
         open(UNPAIR,"> $trim_output_unpaired");
     }
  }

  while ($h1 = <IN>) {  # read first header
        $drop_1=0;
        $raw_seq_num++;
        $s = <IN>;  # read sequence
        chomp $s;
        $h2 = <IN>;  # read second header
        $q = <IN>;  # read quality scores
        chomp $q;
        $len = length($q);
        $total_raw_seq_len += $len;
        $drop_1=1 if ($len < $opt_min_L);
        if ($drop_1==0){
           ($s1_nontrim,$q1_nontrim,$s_trimmed,$q_trimmed,$pos5,$pos3)= &adept_trim($len,$s,$q, $basequality);
           $trim_len=length($s_trimmed);
           $drop_1=1 if ($trim_len < $opt_min_L);
        }
        if ($drop_1==0)
        {
            ($seq_r,$drop_1)=&get_base_and_quality_info($s_trimmed,$q_trimmed,$trim_len,$pos5,$pos3,\%stats);
            %stats=%{$seq_r};
        }
        if ($drop_1==0){  # pass all filters...
            $q_trimmed=&quality_encoding_coversion($q_trimmed,$ascii,$out_offset) if ($ascii != $out_offset);
            $trim_seq_num_1++;
            push @trim_seq_len, $trim_len;
            $total_trim_seq_len += $trim_len;
        }
        
        if ($input2) {
           $drop_2=0;
           $r2_h1 = <IN2>;
           $raw_seq_num++;
           $r2_s = <IN2>;  # mate read sequence
           chomp $r2_s;
           $r2_h2 = <IN2>;  # mate read second header
           $r2_q = <IN2>;  # mate read quality scores
           chomp $r2_q;
           $r2_len = length($r2_q);
           $total_raw_seq_len += $r2_len;
           $drop_2=1 if ($r2_len < $opt_min_L);
           if ($drop_2==0){
              ($s2_nontrim,$q2_nontrim,$r2_s_trimmed, $r2_q_trimmed, $pos5, $pos3)=&adept_trim($r2_len,$r2_s,$r2_q,$basequality);
               $r2_trim_len=length($r2_s_trimmed);
               $drop_2=1 if ($r2_trim_len < $opt_min_L);
           }
           if ($drop_2==0)
           {
               ($seq_r,$drop_2)=&get_base_and_quality_info($r2_s_trimmed,$r2_q_trimmed,$r2_trim_len,$pos5,$pos3,\%stats);
               %stats=%{$seq_r};
           }
          
           if ($drop_2==0){ # pass all filters
               $r2_q_trimmed=&quality_encoding_coversion($r2_q_trimmed,$ascii,$out_offset) if ($ascii != $out_offset);
               $trim_seq_num_2++;
               push @trim_seq_len, $r2_trim_len;
               $total_trim_seq_len += $r2_trim_len;
           }
        }
        if ($drop_1==0 and $drop_2==0)
        {
            $paired_seq_num +=2;
            $total_paired_bases += $trim_len + $r2_trim_len;  
        }
 
        # output trimmed files
        if (! $qc_only)
        {

        if ($nontrimed_output) 
          {
                print OUTQC $h1.$s1_nontrim."\n".$h2.$q1_nontrim."\n";
              if ($input2) {   print OUTQC2 $r2_h1.$s2_nontrim."\n".$r2_h2.$q2_nontrim."\n";}

          }

            if ($drop_1==0 and $drop_2==0)
            {
	        print OUT $h1.$s_trimmed."\n".$h2.$q_trimmed."\n";
                print OUT2 $r2_h1.$r2_s_trimmed."\n".$r2_h2.$r2_q_trimmed."\n";
            }
            elsif ($drop_1==1 and $drop_2==0)
            {
                print DISCARD $h1.$s."\n".$h2.$q."\n" if ($output_discard);
                print UNPAIR $r2_h1.$r2_s_trimmed."\n".$r2_h2.$r2_q_trimmed."\n";
            }
            elsif ($drop_1==0 and $drop_2==1)
            {
                if ($input2)
                {
                   print UNPAIR $h1.$s_trimmed."\n".$h2.$q_trimmed."\n";
                   print DISCARD $r2_h1.$r2_s."\n".$r2_h2.$r2_q."\n" if ($output_discard);
                }
                else
                {
                   print OUT $h1.$s_trimmed."\n".$h2.$q_trimmed."\n";
                }
            }
            elsif ($output_discard)
            {   
                print DISCARD $h1.$s."\n".$h2.$q."\n";
                if ($input2)
                {
                   print DISCARD $r2_h1.$r2_s."\n".$r2_h2.$r2_q."\n";
                }
            }
        }
  } # end while

  my ($trim_seq_len_std,$trim_seq_len_avg,$max,$min,$median)=&standard_deviation(@trim_seq_len);
  
  close(IN);
  close(OUT);
  close IN2 if ($input2);
  close OUT2 if ($input2);
  close UNPAIR if ($input2);
  close DISCARD;
  $stats{raw_seq_num}=$raw_seq_num;
  $stats{trim_seq_num}=$trim_seq_num_1+$trim_seq_num_2;
  $stats{trim_seq_num_1}=$trim_seq_num_1;
  $stats{trim_seq_num_2}=$trim_seq_num_2;
  $stats{total_raw_seq_len}=$total_raw_seq_len;
  $stats{total_trim_seq_len}=$total_trim_seq_len;
  $stats{trim_seq_len_std}=$trim_seq_len_std;
  $stats{trim_seq_len_avg}=$trim_seq_len_avg;
  $stats{max}=$max;
  $stats{min}=$min;
  $stats{median}=$median;
  $stats{paired_seq_num}=$paired_seq_num;
  $stats{total_paired_bases}=$total_paired_bases;
  $stats{trim_file_name_1}=$trim_output_1;
  $stats{trim_file_name_2}=$trim_output_2;
  $stats{numOfReadsWithAdapter}=$numOfReadsWithAdapter;
  return \%stats;
}

     
sub adept_trim
{
    # bwa trim implementation from both 5' and 3' end
    # at least scan 5 bases from both end and 2 more bases after the negative area
    my ($len,$trim_seq,$trim_qual_seq,$basequal) = @_;
    my @trim_seq=split //, $trim_seq;
   my @trim_qual_seq=split //, $trim_qual_seq;
   my %sortbasequal=%{$basequal};
  my @numqual;
  for my $pos(0..$#trim_qual_seq)
  {
   push @numqual, (ord($trim_qual_seq[$pos]) - $ascii);
  }

 
    my $at_least_scan=5;
    my $num_after_neg=2;
    $at_least_scan = $len if $len <5;
    $num_after_neg = $len if $len <2;
    # trim 3' end
    my $pos_3 = $len;
    my $final_pos_3 =$pos_3+1;  
    my $area=0;  
    my $maxArea=0;
    while ($at_least_scan) 
    {
        $at_least_scan--;    
        if ($pos_3>$num_after_neg and $area>=0) {$at_least_scan=$num_after_neg;}
    #	    $area += $opt_q - (ord(substr($q,$pos_3-1,1))-$ascii);
            $area += $opt_q - $numqual[$pos_3-1];
	    if ($area > $maxArea) {
		    $maxArea = $area;
		    $final_pos_3 = $pos_3;
	    }
	    $pos_3--;
    }
    # trim 5' end
    my $pos_5=1; 
    my $final_pos_5=0;
    $maxArea=0;
    $area=0;
    $at_least_scan=5;
    $at_least_scan = $len if $len <5;
    while ($at_least_scan) 
    {
        $at_least_scan--; 
        if ($pos_5<($final_pos_3-$num_after_neg) and $area>=0) {$at_least_scan=$num_after_neg;}
#	    $area += $opt_q - (ord(substr($q,$pos_5-1,1))-$ascii);
            $area += $opt_q - $numqual[$pos_5-1];
	    if ($area > $maxArea) {
		    $maxArea = $area;
		    $final_pos_5 = $pos_5;
	    }
	    $pos_5++;
    }
         foreach  (my $tmpi=0;$tmpi<$final_pos_5;$tmpi++) {$trim_seq[$tmpi]='n';}
         foreach  (my $tmpi=$final_pos_3;$tmpi<$len;$tmpi++) {$trim_seq[$tmpi]='n';}
  
  #ADEPT trimming;
 
  my $tmpcount=0;
  my $tmpstart=0;
  my $longestcount=$minilength;
  my ($longestseq, $longestqual);
  my $prequal2=$numqual[0];
  my $nextqual2=$numqual[-1];
   my $prequal1=$numqual[1];
  my $nextqual1=$numqual[-2];

  my $preavequal2=$sortbasequal{0}[int($cutlimit2*$#{$sortbasequal{0}}+0.5)];
   my $preavequal1=$sortbasequal{1}[int($cutlimit1*$#{$sortbasequal{1}}+0.5)];
  my $nextavequal2=$sortbasequal{$#numqual}[int($cutlimit2*$#{$sortbasequal{$#numqual}}+0.5)];
  my $nextavequal1=$sortbasequal{$#numqual-1}[int($cutlimit1*$#{$sortbasequal{$#numqual-1}}+0.5)];

 my ($uppos, $lowpos);
  for my $pos(0..$#numqual)
  {
    if ($pos>1) {
       $prequal1=$numqual[$pos-1];
     $prequal2=$numqual[$pos-2];
     $preavequal1=$sortbasequal{$pos-1}[int($cutlimit1*$#{$sortbasequal{$pos-1}}+0.5)];
     $preavequal2=$sortbasequal{$pos-2}[int($cutlimit2*$#{$sortbasequal{$pos-2}}+0.5)];
   }
    if ($pos<($#numqual-1)) {

     $nextqual1=$numqual[$pos+1];
      $nextqual2=$numqual[$pos+2];
  $nextavequal1=$sortbasequal{$pos+1}[int($cutlimit1*$#{$sortbasequal{$pos+1}}+0.5)];
      $nextavequal2=$sortbasequal{$pos+2}[int($cutlimit2*$#{$sortbasequal{$pos+2}}+0.5)];
        }
      my $qual =  $numqual[$pos];
        $uppos=int($upperc*$#{$sortbasequal{$pos}}+0.5);
       $lowpos=int($lowperc*$#{$sortbasequal{$pos}}+0.5);

    if ($qual<=$sortbasequal{$pos}[$uppos]  || $trim_seq[$pos] eq 'N'|| $trim_seq[$pos] eq 'n') {

     if ($qual<=$sortbasequal{$pos}[$lowpos]  || $trim_seq[$pos] eq 'N' || $trim_seq[$pos] eq 'n' ) {
        if ($tmpcount>=$minilength) {
     my $tmpseq=substr($trim_seq, $tmpstart, $tmpcount);
     my $tmpqual_seq=substr($trim_qual_seq, $tmpstart, $tmpcount);

  if ($tmpcount>=$longestcount) {
       $longestcount=$tmpcount;
       $longestseq=$tmpseq;
       $longestqual=$tmpqual_seq;
       $final_pos_5=$tmpstart;
     # $final_pos_3=$tmpstart+$tmpcount;
      }
            }
         $trim_seq[$pos] = 'n';
        $tmpstart=$pos+1; $tmpcount=0;
         } else {
  if ( ($qual < $prequal2*$discperc || $qual < $nextqual2*$discperc) && ($qual < $prequal1*$discperc || $qual < $nextqual1*$discperc)
         &&  ($nextqual2 < $nextavequal2 || $prequal2 < $preavequal2) && ($nextqual1 < $nextavequal1 || $prequal1 < $preavequal1) )
   {
                 if ($tmpcount>=$minilength) {
     my $tmpseq=substr($trim_seq, $tmpstart, $tmpcount);
     my $tmpqual_seq=substr($trim_qual_seq, $tmpstart, $tmpcount);

     if ($tmpcount>=$longestcount) {
       $longestcount=$tmpcount;
       $longestseq=$tmpseq;
       $longestqual=$tmpqual_seq;
       $final_pos_5=$tmpstart;
    #  $final_pos_3=$tmpstart+$tmpcount;
      }
            }
           $trim_seq[$pos] = 'n';
        $tmpstart=$pos+1; $tmpcount=0;
          } else {
          $tmpcount++;
          }
        }
     } else { $tmpcount++;}
  }
      my $nontrimseq=join "", @trim_seq;

  #end ADEPT trim 
  if ($longestseq && $longestqual) 
   {
      $final_pos_3=$final_pos_5+length($longestseq);
   return ($nontrimseq,$trim_qual_seq,$longestseq,$longestqual,$final_pos_5,$final_pos_3);
   } else {
       return ($nontrimseq,$trim_qual_seq,'a',0,0,0);
   }
} 

sub get_base_and_quality_info
{
     my $s=shift;
     my $q=shift;
     my $len=shift;
     my $start_pos=shift;
     my $end_pos=shift;
     my $stats=shift;
     my $total_q=0;
     my $avg_q=0;
     my $drop=0;
     my ($a_Base,$t_Base,$c_Base,$g_Base,$n_Base)=(0,0,0,0,0);
     my ($a_Base_percent,$t_Base_percent,$c_Base_percent,$g_Base_percent,$n_Base_percent);
     my $base_percent_string;
     my ($previous_base,$dinucleotide,%di_nucleotide);
     my %seq=%{$stats};
     $start_pos++;
     $end_pos--;
    # my $tmpln=length($s);
     for my $pos($start_pos..$end_pos)
     {
         #print "$q\n$tmpln,$pos,$start_pos,$end_pos,$len\n";
         my $q_digit=ord(substr($q,$pos-$start_pos,1))-$ascii;
         $seq{qual}->{$pos}->{$q_digit}++;
         $total_q += $q_digit; 
         my $base=uc(substr($s,$pos-$start_pos,1));
     #    print "$s\n$pos,$start_pos,$end_pos,'seq'\n";
         $seq{base}->{$pos}->{$base}++;
         $a_Base++ if ($base =~ /A/);  
         $t_Base++ if ($base =~ /T/);  
         $c_Base++ if ($base =~ /C/);  
         $g_Base++ if ($base =~ /G/);  
         $n_Base++ if ($base =~ /N/);  
         if ($previous_base and ($previous_base ne $base))
         {
            $dinucleotide = $previous_base.$base;
            $di_nucleotide{$dinucleotide}++;
         }
         $previous_base=$base;
     }
     $avg_q = int ($total_q/$len);
     if ($drop == 1)  #substract the position matrix by one because the dropped read
     {
         for my $pos($start_pos..$end_pos)
         {
            my $q_digit=ord(substr($q,$pos-$start_pos,1))-$ascii;
            $seq{qual}->{$pos}->{$q_digit}--;
            my $base=uc(substr($s,$pos-$start_pos,1));
            $seq{base}->{$pos}->{$base}--;
         }
         return (\%seq,$drop);
     }
     else
     {
         $a_Base_percent = sprintf("%.2f",$a_Base/$len*100);
         $t_Base_percent = sprintf("%.2f",$t_Base/$len*100);
         $c_Base_percent = sprintf("%.2f",$c_Base/$len*100);
         $g_Base_percent = sprintf("%.2f",$g_Base/$len*100);
         $n_Base_percent = sprintf("%.2f",$n_Base/$len*100);
         $seq{ReadAvgQ}->{$avg_q}->{readsNum}++;
         $seq{ReadAvgQ}->{$avg_q}->{basesNum}+=$len;
         $seq{ReadLen}->{$len}++;
         $seq{Base_content}->{A}->{$a_Base_percent}++;
         $seq{Base_content}->{T}->{$t_Base_percent}++;
         $seq{Base_content}->{C}->{$c_Base_percent}++;
         $seq{Base_content}->{G}->{$g_Base_percent}++;
         $seq{Base_content}->{N}->{$n_Base_percent}++;
         $seq{Base_content}->{GC}->{$c_Base_percent+$g_Base_percent}++;
 
         return (\%seq,$drop);
     }
} 

sub standard_deviation {
  my(@numbers) = @_;
  #Prevent division by 0 error in case you get junk data
  return undef unless(scalar(@numbers));
  @numbers= sort {$b<=>$a} @numbers;
  my $max=$numbers[0];
  my $min=$numbers[-1];
  my $median;
  if (scalar(@numbers) % 2) {
      $median = $numbers[int(scalar(@numbers)/2)];
  } else {
      $median = ($numbers[scalar(@numbers)/2] + $numbers[scalar(@numbers)/2 - 1]) / 2;
  }

  # Step 1, find the mean of the numbers
  my $total1 = 0;
  foreach my $num (@numbers) {
  $total1 += $num;
  }
  my $mean1 = $total1 / (scalar @numbers);

  # Step 2, find the mean of the squares of the differences
  # between each number and the mean
  my $total2 = 0;
  foreach my $num (@numbers) {
  $total2 += ($mean1-$num)**2;
  }
  my $mean2 = $total2 / (scalar @numbers);

  # Step 3, standard deviation is the square root of the
  # above mean
  my $std_dev = sqrt($mean2);
  return ($std_dev,$mean1,$max,$min,$median);
}

sub random_subsample
{
  my $total_num=shift;
  my $output_num=shift;
  my %get;
  if ($output_num <= 0)
  {
     $get{1}=0;
  }
  else
  { 
    while (1)
    {
      $get{int(rand($total_num))}=1;
      my $num_of_random=scalar (keys %get);
      last if ($num_of_random == $output_num or $num_of_random == $total_num);
    }
  }
  return (\%get);
}
sub quality_encoding_coversion 
{  
    # given quality acsii string, input offset and output offset
    my $q_string=shift;
    my $input_offset=shift;
    my $out_offset=shift;
    $q_string=~ s/(\S)/chr(ord($1)-$input_offset+$out_offset)/eg;
    return($q_string);
}


sub checkQualityFormat {
    # $offset_value=&checkQualityFormat($fastq_file)
    # $offset_value = -1 means not proper fastq format.
    my $fastq_file=shift;
    # open the files
    if ($fastq_file =~ /gz$/)
    {
        open FQ, "zcat $fastq_file | " or die $!;
    }
    else
    {
        open FQ, "<", $fastq_file or die $!;
    }

    # initiate
    my @line;
    my $l;
    my $number;
    my $offset;
    # go thorugh the file
    my $first_line=<FQ>;
    if ($first_line !~ /^@/) {$offset=-1; return $offset;}
    OUTER:while(<FQ>){
      
      # if it is the line before the quality line
      if($_ =~ /^\+/){
    
       $l = <FQ>; # get the quality line
       @line = split(//,$l); # divide in chars
       for(my $i = 0; $i <= $#line; $i++){ # for each char
        $number = ord($line[$i]); # get the number represented by the ascii char
      
        # check if it is sanger or illumina/solexa, based on the ASCII image at http://en.wikipedia.org/wiki/FASTQ_format#Encoding
        if($number > 74){ # if solexa/illumina
          $offset=64;
          #die "This file is solexa/illumina format\n"; # print result to terminal and die
          last OUTER; 
        }elsif($number < 59){ # if sanger
          $offset=33;
          #die "This file is sanger format\n"; # print result to terminal and die
          last OUTER;
        }
       }
      }
    }
   return $offset;
}

sub split_fastq {
   my ($input,$outDir,$subfile_size)=@_;

   my $seqThisFile = 0;
   my $fileCount = 0;
   my $name="";
   my $seq="";
   my $q_name="";
   my $qual_seq="";
   my @subfiles;
   my $seq_num=0;
   my $total_seq_length=0;
   my ($file_name, $i_path, $i_suffix_1)=fileparse("$input", qr/\.[^.]*/);
   ($file_name, $i_path, $i_suffix)=fileparse("$file_name", qr/\.[^.]*/) if ($i_suffix_1 =~ /gz$/);
   my $file_ext = sprintf("%05d", $fileCount);
   my $out_file_name= $outDir . "/". $file_name . "_" . $file_ext . ".fastq";
   push @subfiles, $out_file_name;
   open (OUTFILE, ">" . $out_file_name ) or die ("Cannot open file for output: $!");
   if ($i_suffix_1 =~ /gz/)
   {
      open (IN, "zcat $input | ") || die "cannot open file $input:$!";
   }
   else
   {
      open (IN, $input) || die "cannot open file $input:$!";
   }
   while (<IN>)
   {
          $name = $_;
          $seq=<IN>;
          $seq =~ s/\n//g;
          while ($seq !~ /\+/)
          {
             $seq .= <IN>;
             $seq =~ s/\n//g;
          }
          my $q_name_pos=index($seq,"+");
          $q_name = substr($seq,$q_name_pos);
          $seq = substr($seq, 0, $q_name_pos);
          my $seq_len = length $seq;
          $qual_seq=<IN>;
          $qual_seq =~ s/\n//g;
          my $qual_seq_len = length $qual_seq;
          while ( $qual_seq_len < $seq_len )
          {
              last if ( $qual_seq_len == $seq_len);
              $qual_seq .= <IN>;
              $qual_seq =~ s/\n//g;
              $qual_seq_len = length $qual_seq;
          }
          $total_seq_length += $seq_len;
          print (OUTFILE "$name$seq\n$q_name\n$qual_seq\n");
          $seq_num++;
          $seqThisFile++;

          if ($seqThisFile == $subfile_size)
          {
              $fileCount++;
	      $seqThisFile = 0;
	      close (OUTFILE) or die( "Cannot close file : $!");
	      $file_ext = sprintf("%05d", $fileCount);
	      $out_file_name= $outDir . "/". $file_name . "_" . $file_ext . ".fastq";
              open (OUTFILE, ">" .  $out_file_name ) or die ("Cannot open file for output: $!") if (! eof);
              push @subfiles, $out_file_name if (! eof);
          }

   }
   my $average_len = $total_seq_length/$seq_num;
   if ( $average_len < $opt_min_L) { print "The input average length $average_len < minimum cutoff length(opt_min_L) $opt_min_L\n."; exit;}
   close (IN)  or die( "Cannot close file : $!");
   close (OUTFILE) or die( "Cannot close file : $!") if (! eof OUTFILE);
   return ($seq_num,@subfiles);
}
   

sub ReverseComplement{
        my $dna = $_[0];
        my $ReverseCompSeq = reverse ($dna);
        $ReverseCompSeq =~ tr/atgcrywsmkATGCRYWSMK/tacgyrswkmTACGYRSWKM/;
        return($ReverseCompSeq);
}



1;


