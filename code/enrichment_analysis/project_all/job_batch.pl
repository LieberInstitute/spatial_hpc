use warnings;
use Getopt::Long;

GetOptions ("jfile=s" => \$infile,
		 "nrun=i"   => \$nrun,   
              "n=i"   => \$n)
or die("Error in command line arguments\n");

$i=0;

open(IN, $infile);
while(<IN>){
	chomp;
	$i++;
	$start=($n-1)*$nrun + 1;
	$end=$n*$nrun;
	if($i >= $start and $i <= $end){
		$command=$_;
		system($command);
	}
	if($i > $end){
		last;
	}
}
close(IN);
	