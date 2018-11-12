### 1 ### Setup automatic Lisa access in Mac terminal
### 2 ### Setup individiual python environment in Lisa
### 3 ### Singularity commands for preprocessing on Lisa

### 1 ###
### HELPER: Terminal mapping to lisa
alias lisa_home='sshfs -p 22 mfailing@lisa.surfsara.nl:/home/ ~/surfsara/ -oauto_cache,reconnect,defer_permissions,noappledouble,negative_vncache'

### 2 ###
# Need specific python version?
module load python # default
module load python/2.7.9
module load python/3.4.2 # used through python3 exectuable

# Need specific package?
module load python # load python
pip install --user nistats # install specific package

### 3 ###

## We can simply copy&paste the following after starting up python

import os

subjects = [1]

batch_string = """# shell for the job:
#PBS -S /bin/bash
#PBS -lwalltime=100:00:00 -lnodes=1:mem64gb
# job requires at most 10 hours, 0 minutes
#     and 0 seconds wallclock time
# call the programs
echo "Job $PBS_JOBID started at `date`" | mail $USER -s "Job $PBS_JOBID"

PYTHONPATH="" singularity run /home/mfailing/poldracklab_fmriprep_latest-2018-06-07-3b37987ebff2.img \
                /home/mfailing/fMRI_NRosT/BIDS /home/mfailing/fMRI_NRosT/BIDS/derivatives participant \
                --participant-label sub-SJ_NR --output-space T1w template fsaverage6 fsaverage5 --nthreads 15 --omp-nthreads 15 --use-syn-sdc --low-mem \
                --fs-license-file /home/mfailing/bin/freesurfer/license.txt -w /scratch

wait          # wait until programs are finished

echo "Job $PBS_JOBID finished at `date`" | mail $USER -s "Job $PBS_JOBID"
"""

basedir = '/home/mfailing/'
os.chdir(basedir)

for subject in subjects:
    working_string = batch_string.replace('SJ_NR', str(subject).zfill(2))
    js_name = os.path.join(basedir, str(subject).zfill(2) + '_nofs.sh')
    of = open(js_name, 'w')
    of.write(working_string)
    of.close()
    print('submitting ' + js_name + ' to queue')
    print(working_string)
    os.system('qsub ' + js_name)




    # Check job queue
    showq -u mfailing
    # Inspect job
    pbs_jobmonitor [jobid] [nodeid]