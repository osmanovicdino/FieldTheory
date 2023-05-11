#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=12:00:00,h_data=2G
## Modify the parallel environment
## and the number of cores as needed:
#$ -t 1-84:1

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
# module load gcc/7.5.0
# module load fftw3/3.3.8-gcc
module load gcc/10.2.0
module load fftw/3.3.9


## substitute the command to run your code
## in the two lines below:
##echo '/usr/bin/time -v hostname'
##/usr/bin/time -v hostname
filename=~/FieldTheory/Code/Hoffman/paramsI.dat
basedir="Invasion3D"
if [ -e ${filename}   ]; then
   # use the unix command sed -n ${line_number}p to read by line
   c0=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $1}'`
   c1=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $2}'`
   l=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $3}'`
   p0=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $1}'`
   A0=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $2}'` 
   dirc=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $3}'`
   fil1=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $4}'`
   echo "read file correctly" 
else
   p0=0.5;
   A0=0.5;
   echo "did not read file correctly"
fi
dirwemake="paramset=${SGE_TASK_ID}"
mkdir /u/scratch/d/dinoo/${basedir}/${dirwemake}
cp ~/FieldTheory/Code/mainInvasion3DImport.cpp /u/scratch/d/dinoo/${basedir}/${dirwemake}
g++ ~/FieldTheory/Code/mainInvasion3DImport.cpp -lm -lfftw3 -L/usr/local/lib/lfftw3.a -std=c++17 -o /u/scratch/d/dinoo/${basedir}/${dirwemake}/angron
cd /u/scratch/d/dinoo/${basedir}/${dirwemake}
./angron $c0 $c1 $l $p0 $A0 $dirc $fil1 > log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####