use warnings;

# change path below to your correct path
$ldsc="/dcs04/lieber/shared/statsgen/LDSC/base/scripts/ldsc.py";
$referenceDir="/dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_EUR_Phase3_plink/";
$hapmap3="/dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/hapmap3_snps/";

open(OUT, ">ldsc_score_jobs.txt");

opendir($dh, "./");
@outdir = grep { !/^\./ && /^out_/ && -d "$_"} readdir($dh);
closedir $dh;

foreach $out (@outdir){
  foreach $i(1..22){
	$bfile=$referenceDir."1000G.EUR.QC.".$i;
	$anno="./$out/chr.".$i.".annot.gz";
	$hapmapSNP=$hapmap3."hm.".$i.".snp";
	print OUT "python $ldsc --l2 --thin-annot --ld-wind-cm 1 ";
	print OUT "--bfile $bfile ";
	print OUT "--anno $anno ";
	print OUT "--out ./$out/chr.$i ";
	print OUT "--print-snps $hapmapSNP\n";
  }
}

close(OUT);