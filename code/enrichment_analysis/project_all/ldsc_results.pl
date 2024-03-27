use warnings;
use Cwd;

$cwd = cwd();

# trait meta
open(IN, "/dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain/trait_names_keys");
while(<IN>){
	chomp;
	@tokens=split(/\t/,$_);
	$filename=$tokens[0];
	$filename =~ s/\.gz/\.out\.results/;
	#$type{$filename}=$tokens[1];
	$trait{$filename}=$tokens[2];
	print "$filename\t$trait{$filename}\n";
}
close(IN);

opendir($dh, "./");
@outdirs = grep { !/^\./ && /_2$/ && -d "$_"} readdir($dh);
closedir $dh;

$count=0;

open(OUT, ">ldsc_results.txt");
print OUT "cell\ttrait\tProp._SNPs	Prop._h2	Prop._h2_std_error	Enrichment	Enrichment_std_error	Enrichment_p	Coefficient	Coefficient_std_error	Coefficient_z-score\n";
foreach $dir (@outdirs){
	$count++;
	print "$count\n";
	$group=$dir;
	$group =~ s/out_//;
	$group =~ s/_2$//;
	chdir $dir;
	foreach $file(<*results>){
		if(not defined $trait{$file}){
			print "$file\n";
		}
		print OUT "$group\t$trait{$file}\t";
		open(IN,$file);
		<IN>;
		$line=<IN>;
		chomp $line;
		@tokens=split(/\t/,$line);
		shift @tokens;
		print OUT join("\t", @tokens), "\n";
		close(IN);
	}
	chdir "..";
}
close(OUT);