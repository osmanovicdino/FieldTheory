#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=3:00:00,h_data=2G
## Modify the parallel environment
## and the number of cores as needed:
#$ -t 1-38:1

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
filename=~/FieldTheory/Code/Hoffman/params7.dat
mydir="TwoFieldTheory5"
if [ -e ${filename}   ]; then
   # use the unix command sed -n ${line_number}p to read by line
    x12=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $1}'` 
    x13=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $2}'` 
    x23=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $3}'` 
    den1=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $4}'` 
    den2=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $5}'` 
    cut=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $6}'` 
   echo "read file correctly" 
else
   x12=1.0;
   x13=1.0;
   x23=1.0;
   den1=0.0;
   den2=0.0;
   cut=100;
   echo "did not read file correctly"
fi
dirwemake="paramset${SGE_TASK_ID}"
mkdir /u/scratch/d/dinoo/${mydir}/${dirwemake}
cp ~/FieldTheory/Code/mainCHDoubleInject.cpp /u/scratch/d/dinoo/${mydir}/${dirwemake}
g++ ~/FieldTheory/Code/mainCHDoubleInject.cpp -lm -lfftw3 -L/usr/local/lib/lfftw3.a -std=c++17 -o /u/scratch/d/dinoo/${mydir}/${dirwemake}/angron
cd /u/scratch/d/dinoo/${mydir}/${dirwemake}
./angron $x12 $x13 $x23 $den1 $den2 $cut > log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####