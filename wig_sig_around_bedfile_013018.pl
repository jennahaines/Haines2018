#########################################################################################
##########   wig_sig_around_bedfile_013018.pl                                 ###########
##########   summarizes wig signal around a bedfile                           ###########
##########   Input - bed file, wig file directory                             ###########
##########   Output - summarized file                                         ###########
#########################################################################################
#!/usr/bin/perl

use warnings;

#average signal over a bed region
my $peaks =shift; #bed file 
my $dir = shift;  #wig file dir
# my $NF = shift;  #normalization table
my $output = shift;
#my $win = shift;

# my %NF_t;

# open (NFT,'<', $NF);
# while (<NFT>){
#  next if(/#/);

#  my ($FS, $v) = split /\s+/,$_;
# $NF_t{$FS} = $v;

# }

my %HOAg=();

open (OUTPUT, '>', $output);
#print OUTPUT "#gen_ann,average signal within cluster region\n";
#print OUTPUT "#ID\tchr\tpeak_pos1\tpeak_pos2\tinterval_length\t";  #for union of peaks file
open (INPUTp, '<', $peaks);

my $N=1;
while (<INPUTp>){
   if (/^chr(\d|X)/) {
    @split = split (/\s+/, $_);

   # my $PID = $split[0];
    $HOAg{$N} = [@split];
    $N +=1;
  } else {
     @split = split (/\s+/, $_);
     $st = join ("\t",@split);
     print OUTPUT "$st";
    }
}


opendir DIR, $dir or die "can't open $dir: $!";

while ( my $file = readdir DIR ) {
  if ($file =~ /(.+).wig$/) {
    $HM = "$1";
    print OUTPUT "\t$HM";

    $In_file = "$dir/$file";
    print "infile:$file\n";

    my %HOA;
    my $chr2 = "NA";
    my %chr_sz=(
	   "chr2L" => 23011344, #new size info, aug 2014
	   "chr2R" => 21146608,
	   "chr3L" => 24543457,
	   "chr3R" => 27905053,
	   "chr4" => 1351757,
	   "chrX" => 22332727);

    foreach $chr (keys %chr_sz){
      $n = int($chr_sz{$chr}/10);
      @a = (0)x$n;
      $HOA{$chr} = [@a];
		}


    open (INPUT, '<', $In_file);
    while ( <INPUT> ) {
      if ($_ =~/^tr/){
	next}elsif ($_ =~ /chrom=(\w+)\s+/){
	  $chr2 = $1;
	}elsif ($_ =~/^\d+/){
	  @s=split (/\s+/, $_);
	  $p =int( $s[0]/10);		      		     
	  $HOA{$chr2}[$p] += $s[1];		     
	}
}
close INPUT;

    my $total=0;
    my $L;
    foreach $i (keys %HOAg){
      $chr3 = $HOAg{$i}[0];  #
  
   $st_p = int($HOAg{$i}[1]/10); # 5873439 becomes 587343
   $en_p = int($HOAg{$i}[2]/10);

    $L = ($en_p - $st_p);


      @ar= @{$HOA{$chr3}}[$st_p..$en_p];
      #print "@ar\n";

      $total += $_ for @ar;
  $val = $total/$L;
  #$vf= sprintf("%.2f", $val);
  push @{$HOAg{$i}}, $val;
   #push @{$HOAg{$i}}, $total;
	  
      $total=0
    }
  }
}

print OUTPUT "\n";

foreach $i (sort {$a <=> $b}(keys %HOAg)){
  $str = join ("\t", @{$HOAg{$i}});
print OUTPUT "$str\n"
}
close DIR;
close OUTPUT;