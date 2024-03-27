# simpleArray.pl, created by Shizhong Han, Oct-20-2023
# This is a script that combines array tasks with a
# perl script to process many short runs. Array jobs are convenient
# for running lots of tasks, but if each task is short, they
# quickly become inefficient, taking more time to schedule than
# they spend doing any work and bogging down the scheduler for
# all users. this script applies to shared partition in jhpce3

# usage: perl simpleArray.pl <--job jobfile> [--nrun #] <--out batchfile>

# arguments:
# -j, or --job	input job file
# -n, or --nrun	the number of short runs that each SLURM task should do
# -o, or --out	output batch script

# example: 
# simpleArray.pl -j plink_job.txt -n 40 -o plink_jobs.batch

use warnings;
use Getopt::Long;

GetOptions(
    'job=s'       => \$job,
    'nrun=i'       => \$nrun,
    'out=s'       => \$out    
) or die "Incorrect usage!\n";

# perl batch script
$job_batch="job_batch.pl";

# directory for err_out
mkdir $oe=$job."_err_out";

# job files
open (IN, $job);
$i=0;
while(<IN>){
	$i++;	
}
close(IN);

# number of SLURM tasks
$ntasks=$i/$nrun;
if (not $ntasks =~ /^-?\d+$/){
	$ntasks = int ($i/$nrun) + 1;
}

# batch file
open (OUT, ">$out");

print OUT "#!/bin/bash\n";
print OUT "#SBATCH --job-name=$job\n";
print OUT "#SBATCH --output=./$oe/jobs_%A_%a.out\n";
print OUT "#SBATCH --error=./$oe/jobs_%A_%a.err\n";
print OUT "#SBATCH --partition=shared\n";
print OUT "#SBATCH --mem-per-cpu=5G --cpus-per-task=1\n";
print OUT "#SBATCH --array=1-${ntasks}%${ntasks}\n\n\n";

print OUT "echo \"**** Job starts ****\"\n";
print OUT "date\n";
print OUT "module load conda_R/4.0\n";

print OUT "perl $job_batch --jfile $job --nrun $nrun --n \$SLURM_ARRAY_TASK_ID\n";

print OUT "echo \"**** Job ends ****\"\n";
print OUT "date\n";

close(OUT);





