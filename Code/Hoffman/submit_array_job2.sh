#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=6:00:00,h_data=2G
## Modify the parallel environment
## and the number of cores as needed:
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea
#$ -t 1-7:1

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load gcc/7.5.0
module load fftw3/3.3.8-gcc

## substitute the command to run your code
## in the two lines below:
##echo '/usr/bin/time -v hostname'
##/usr/bin/time -v hostname
filename=~/FieldTheory/Code/Hoffman/params2.dat
mydir="FieldTheory2"
if [ -e ${filename}   ]; then
   # use the unix command sed -n ${line_number}p to read by line
   D0=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $1}'` 
   echo "read file correctly" 
else
   D0=1.0;
   echo "did not read file correctly"
fi
dirwemake="D0=${D0}"
mkdir /u/scratch/d/dinoo/${mydir}/${dirwemake}
cp ~/FieldTheory/Code/maintest4.cpp /u/scratch/d/dinoo/${mydir}/${dirwemake}
g++ ~/FieldTheory/Code/maintest4.cpp -lm -lfftw3 -L/usr/local/lib/lfftw3.a -std=c++11 -o /u/scratch/d/dinoo/${mydir}/${dirwemake}/angron
cd /u/scratch/d/dinoo/${mydir}/${dirwemake}
./angron $D0 > log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####