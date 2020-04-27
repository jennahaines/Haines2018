##########################################################################################
##########   xl-bed2wig-dirproc-HOA-2014-new-chrsz                              ##########
##########   Input: Bed files                                                   ##########
#########    Output: Wig files                                                  ##########
##########################################################################################


#!/usr/bin/perl

use warnings;

my $dir = shift;

opendir DIR, $dir or die "can't open $dir: $!";
while ( my $file = readdir DIR ) {
  if ($file =~ /^(.+).(BED|bed|bedgraph)/) {
    $In_file = "$dir/$file";
    print "infile:$In_file\n";
    my $output ="$dir/$1.wig";
    print "output:$output\n";
    my $wigname = "$1";
    my %HOA = ();
    my %chr_sz=(
                "chr2L" => 23011344, #new size info, aug 2014
                "chr2R" => 21146608,
                "chr3L" => 24543457,
                "chr3R" => 27905053,
                "chr4" => 1351757,
                "chrX" => 22332727 );

    foreach $chr (keys %chr_sz){
      $n = int($chr_sz{$chr}/10);
      @a = (0)x$n;
      $HOA{$chr} = [@a];
    }

    my ($chr2, $stp, $endp);
    open (INPUT, '<', $In_file);
    while (<INPUT>) {
      next if (/(tr|U|Het|dmel)/);
      my @spl = split /\s+/, $_;
      # $chr2 = $spl[0];
      # $stp = $spl[1]; 
      # $endp = $spl[2];
      my @n;
      my $chp;
      $chr_n = $spl[0];
      $chr_n =~ tr/chr//d;
      $chr2 = "chr$chr_n";
        
      if ($spl[2] < $chr_sz{$chr2}) {
        foreach $n($spl[1] .. $spl[2]) {                      
          $chp = int($n/10);
          $HOA{$chr2}[$chp] +=1;
        }
      }
    }                           
                        
    close INPUT;

    open (OUTPUT, '>', $output);
    print OUTPUT "track type=wiggle_0 name=\"$wigname\" description=\"$wigname\"visibility=full autoScale=off maxHeightPixels=100:50:20\n";
    foreach $key (keys %HOA){
      print OUTPUT "variableStep chrom=$key span=10\n";
      foreach $i (0 .. $#{$HOA{$key}}) {
        my $m = ($i+1)*10;
        my $v = $HOA{$key}[$i];
        if ($v > 0){
          print OUTPUT "$m\t$v\n";
                }
      }
    }  
    close OUTPUT;
  }
}
closedir(DIR);