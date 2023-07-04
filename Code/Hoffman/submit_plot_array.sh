#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=3:00:00,h_data=8G
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
# module load gcc/10.2.0
# module load fftw/3.3.9
module load python/3.9.6
module load ffmpeg

## substitute the command to run your code
## in the two lines below:
##echo '/usr/bin/time -v hostname'
##/usr/bin/time -v hostname

dirwemake="chemistry${SGE_TASK_ID}"
cd /u/scratch/d/dinoo/WavesSingle/${dirwemake}
export OMP_NUM_THREADS=1
num=`ls field1*.csv | wc -l`
if test $num -ge 100;
then
    for file in field0*.csv; do python3 ~/FieldTheory/Code/Plotting/PlotDirectory.py "$file"; done
    ffmpeg -pattern_type glob -i 'field0*.png' -s 1920x1080 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" test.mp4
    for file in field1*.csv; do python3 ~/FieldTheory/Code/Plotting/PlotDirectory.py "$file"; done
    ffmpeg -pattern_type glob -i 'field1*.png' -s 1920x1080 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" test2.mp4
fi
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####