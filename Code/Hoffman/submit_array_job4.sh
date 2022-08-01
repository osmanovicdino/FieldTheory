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
#$ -t 1-125:1

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
filename=~/FieldTheory/Code/Hoffman/params3.dat
mydir="TwoFieldTheory1"
if [ -e ${filename}   ]; then
   # use the unix command sed -n ${line_number}p to read by line
    x12=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $1}'` 
    x13=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $2}'` 
    x23=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $3}'` 
   echo "read file correctly" 
else
   x12=1.0;
   x13=1.0;
   x23=1.0;
   echo "did not read file correctly"
fi
dirwemake="paramset${SGE_TASK_ID}"
mkdir /u/scratch/d/dinoo/${mydir}/${dirwemake}
cp ~/FieldTheory/Code/mainCHDouble.cpp /u/scratch/d/dinoo/${mydir}/${dirwemake}
g++ ~/FieldTheory/Code/mainCHDouble.cpp -lm -lfftw3 -L/usr/local/lib/lfftw3.a -std=c++11 -o /u/scratch/d/dinoo/${mydir}/${dirwemake}/angron
cd /u/scratch/d/dinoo/${mydir}/${dirwemake}
./angron $x12 $x13 $x23 > log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####