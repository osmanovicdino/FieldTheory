#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=9:00:00,h_data=2G
## Modify the parallel environment
## and the number of cores as needed:
#$ -t 1-17:1

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load gcc/10.2.0
module load fftw/3.3.9

## substitute the command to run your code
## in the two lines below:
##echo '/usr/bin/time -v hostname'
##/usr/bin/time -v hostname
filename=~/FieldTheory/Code/Hoffman/allfiles.dat
if [ -e ${filename}   ]; then
   # use the unix command sed -n ${line_number}p to read by line
   wt=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $1}'`
   echo "read file correctly" 
else
   wt=res1.csv
   echo "did not read file correctly"
fi
dirwemake="chemistry${SGE_TASK_ID}"
ftdir="WavesSingle"
mkdir /u/scratch/d/dinoo/${ftdir}/${dirwemake}
cp ~/FieldTheory/Code/mainFractionalDiffusionNormalChem.cpp /u/scratch/d/dinoo/${ftdir}/${dirwemake}
cp ~/FieldTheory/Code/InitialConditions/Cond6/${wt} /u/scratch/d/dinoo/${ftdir}/${dirwemake}/res.csv
g++ ~/FieldTheory/Code/mainFractionalDiffusionNormalChem.cpp -lm -lfftw3 -L/usr/local/lib/lfftw3.a -std=c++17 -o /u/scratch/d/dinoo/${ftdir}/${dirwemake}/angron
cd /u/scratch/d/dinoo/${ftdir}/${dirwemake}
./angron 'res.csv' > log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####