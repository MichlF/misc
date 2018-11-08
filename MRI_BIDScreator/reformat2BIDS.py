'''
Some doc string
@Michel Failing, 2018
'''

import os, re, json, pandas, B0_separation, ees_function
from shutil import copyfile
from pprint import pprint

def folderCreator(path):
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            raise
        pass

### Start script
# FOR NOW:
sesIdx = '01'
nrRuns = 9
no_trialsBlock = 24
time_memEvent = 7.7
time_beginMemEvent = 10.5
time_trial = 19.6

# What to do?
construct = True
behavior = True

# Paths
pathBEHdata = '/Users/michlf/Dropbox/Work/Data/fMRI/NegativeTemplate/forScanner/mri/beh'
pathMRIdata = '/Users/michlf/NROST_working/' 
pathBIDS = '/Users/michlf/NROST_working/fMRI_NROST/' #'/Users/michlf/Documents/GitHub/fMRI_NRoST

# Create directories
if any('sub-' in s for s in os.listdir(pathMRIdata)):  # check whether there is any subject data
    # Load in the config file
    with open('/Users/michlf/NROST_working/config.json') as f:
        data = json.load(f)
    pprint(data)

    # Some other preparations
    os.chdir(pathMRIdata)
    folderlist = os.listdir(pathMRIdata)
    folderlist.sort()
    folderCreator(pathBIDS+'code/')

    # Create folder structure for any subject found
    subList = []
    for i in range(len(folderlist)):
        try:
            subList.append(folderlist[i].split("sub-",1)[1])
            nrSubj = folderlist[i].split("sub-",1)[1]
            print('\n...running subject {0}...\n'.format(nrSubj))
        except:
            continue # if another file, move on to next iteration
        pathSubj = pathBIDS + 'sub-' + nrSubj
        folderCreator(pathSubj+'/ses-{0}/anat/'.format(sesIdx))
        folderCreator(pathSubj+'/ses-{0}/func/'.format(sesIdx))
        folderCreator(pathSubj+'/ses-{0}/fmap/'.format(sesIdx))
        folderCreator(pathSubj+'/ses-{0}/beh/'.format(sesIdx))

        ### Construct name, then copy and rename files
        if construct:
            # Anatomy (T1w, T2)
            path2Original = pathMRIdata + folderlist[i] + '/' + [s for s in os.listdir(pathMRIdata+folderlist[i]) if "_3DT1w" in s and ".nii.gz" in s][0]
            path2New = pathSubj + '/ses-{1}/anat/sub-{0}_ses-{1}_T1w.nii.gz'.format(nrSubj, sesIdx)
            copyfile(path2Original, path2New)
            path2Original = pathMRIdata + folderlist[i] + '/' + [s for s in os.listdir(pathMRIdata+folderlist[i]) if "_3DT2" in s and ".nii.gz" in s][0]
            path2New = pathSubj + '/ses-{1}/anat/sub-{0}_ses-{1}_T2w.nii.gz'.format(nrSubj, sesIdx)
            copyfile(path2Original, path2New)
            
            # Fmap and func from all runs including the localizer
            # B0 separation
            path2Original = pathMRIdata + folderlist[i] + '/' + [s for s in os.listdir(pathMRIdata+folderlist[i]) if "B0" in s and ".nii.gz" in s][0]
            B0_separation.separation(path2Original, pathSubj+'/ses-{0}/fmap/'.format(sesIdx), nrSubj, 1)
            
            # EPI (fmap) and BOLD (func)
            path_intendedList = []
            path_fmapDir = pathSubj + '/ses-{0}/fmap/'.format(sesIdx)
            for runIdx in range(1,nrRuns+1): # 8 task runs + 1 localizer

                # Task
                if runIdx<nrRuns:
                    try:
                        # EPI
                        path2Original = pathMRIdata + folderlist[i] + '/' + [s for s in os.listdir(pathMRIdata+folderlist[i]) if "{0}_epi_".format(runIdx) in s and ".nii.gz" in s][0]
                        path2New = pathSubj + '/ses-{2}/fmap/sub-{0}_ses-{2}_dir-NR_run-{1}_epi.nii.gz'.format(nrSubj, runIdx, sesIdx)
                        copyfile(path2Original, path2New)
                        # BOLD
                        path2Original = pathMRIdata + folderlist[i] + '/' + [s for s in os.listdir(pathMRIdata+folderlist[i]) if "{0}_bold_".format(runIdx) in s and ".nii.gz" in s][0]
                        path2New = pathSubj + '/ses-{2}/func/sub-{0}_ses-{2}_task-NRoST_run-{1}_bold.nii.gz'.format(nrSubj, runIdx, sesIdx)
                        copyfile(path2Original, path2New)
                        # JSON
                        path_fmap_PAR = pathMRIdata + folderlist[i] + '/' + [s for s in os.listdir(pathMRIdata+folderlist[i]) if "{0}_bold_".format(runIdx) in s and ".PAR" in s][0]
                        path_intended = '"func/' + 'sub-{0}_ses-{2}_task-NRoST_run-{1}_bold.nii.gz"'.format(nrSubj, runIdx, sesIdx)
                        path_fmap = path_fmapDir + 'sub-{0}_ses-{2}_dir-NR_run-{1}'.format(nrSubj, runIdx, sesIdx)
                    except:
                        print('runIdx {0} did not work for subject {1}. File missing?'.format(runIdx,nrSubj))

                # Localizer
                elif runIdx==nrRuns:
                    try:
                        # EPI
                        path2Original = pathMRIdata + folderlist[i] + '/' + [s for s in os.listdir(pathMRIdata+folderlist[i]) if "{0}_epi_".format('Local') in s and ".nii.gz" in s][0]
                        path2New = pathSubj + '/ses-{2}/fmap/sub-{0}_ses-{2}_dir-NR_run-{1}_epi.nii.gz'.format(nrSubj, runIdx, sesIdx)
                        copyfile(path2Original, path2New)
                        # BOLD
                        path2Original = pathMRIdata + folderlist[i] + '/' + [s for s in os.listdir(pathMRIdata+folderlist[i]) if "{0}-bold_".format('Local') in s and ".nii.gz" in s][0]
                        path2New = pathSubj + '/ses-{2}/func/sub-{0}_ses-{2}_task-Localizer_run-{1}_bold.nii.gz'.format(nrSubj, runIdx, sesIdx)
                        copyfile(path2Original, path2New)
                        # JSON
                        path_fmap_PAR = pathMRIdata + folderlist[i] + '/' + [s for s in os.listdir(pathMRIdata+folderlist[i]) if "{0}-bold_".format('Local') in s and ".PAR" in s][0]
                        path_intended = '"func/' + 'sub-{0}_ses-{2}_task-Localizer_run-{1}_bold.nii.gz"'.format(nrSubj, runIdx, sesIdx)
                        path_fmap = path_fmapDir + 'sub-{0}_ses-{2}_dir-NR_run-{1}'.format(nrSubj, runIdx, sesIdx)
                    except:
                        print('Localizer run did not work. File missing?')

            ### JSON files
                # Individual EPIs (for each run incl. localizer)
                try:
                    wfs = ees_function.readWfs(path_fmap_PAR)
                    print('water fat shift for run {0}: {1} px'.format(runIdx, wfs))
                    EES = ees_function.calculateEES(path_fmap+'_epi.nii.gz', wfs)
                    File = open(path_fmap + '_epi.json', 'w')
                    File.write(''.join('{')+ '\n')
                    File.write(''.join('    "PhaseEncodingDirection": "j-",')+ '\n')
                    File.write(''.join('    "TotalReadoutTime": ' + '{0}'.format(EES) +',')+ '\n')
                    File.write(''.join('    "IntendedFor": [' + path_intended + ']')+ '\n')
                    File.write(''.join('}')+ '\n')
                    File.close()
                    path_intendedList.append(path_intended) # save for phase difference JSON
                except:
                    print('JSON did not work. File missing?')

            # Phasediff (only needed once, so out of the loop)
            File = open(path_fmapDir + 'sub-{0}_ses-{1}_phasediff.json'.format(nrSubj, sesIdx), 'w')
            File.write(''.join('{')+ '\n')
            File.write(''.join('    "EchoTime1": 0.003,')+ '\n')
            File.write(''.join('    "EchoTime2": 0.008,')+ '\n')
            for run in range(nrRuns):
                if run == 0:
                    File.write(''.join('    "IntendedFor": [' + path_intendedList[run] + ',')+ '\n')
                elif run == nrRuns-1:
                    File.write(''.join('                    ' + path_intendedList[run] + ']')+ '\n')
                else:
                    File.write(''.join('                    ' + path_intendedList[run] + ',')+ '\n')
            File.write(''.join('}')+ '\n')
            File.close()

        ### Physiology (still to do)
        pass

        ### Behavior
        if behavior:
            # Experimental run
            path2Original = pathBEHdata + '/' + [s for s in os.listdir(pathBEHdata) if "subject-{0}_".format(int(nrSubj)) in s and ".csv" in s][0]
            df = pandas.read_csv(path2Original)
            # Load beh and create/change a few variables to comply with BIDS
            df['duration'] = time_memEvent
            df['trial_type'] = 'memory'
            df['responseTime'] = round(df['responseTime'])
            df['condition'] = df['cond_template'].map(str) + '-' + df['cond_category']
            # Now, chop up the file in individual runs
            for runIdx in range(1,nrRuns):
                # Cut the runs
                df_final = df.truncate(before=(runIdx-1)*no_trialsBlock, after=no_trialsBlock-1+((runIdx-1)*no_trialsBlock))
                if int(nrSubj) < 7: # One TR at the end of the first trial of each run is missing
                    df_final['onset'] = [round(num * time_trial + time_beginMemEvent - 0.7, 2) for num in range(no_trialsBlock)] # adjust all onsets by removing 1 TR
                    df_final['onset'][time_beginMemEvent] = 10.5 # change the first onset
                else:
                    df_final['onset'] = [round(num * time_trial + time_beginMemEvent, 2) for num in range(no_trialsBlock)]
                # Save
                path2New = pathSubj + '/ses-{2}/func/sub-{0}_ses-{2}_task-NRoST_run-{1}_events.tsv'.format(nrSubj, runIdx, sesIdx)
                df_final.to_csv(path2New, sep='\t', columns=['onset', 'duration', 'trial_type', 'responseTime', 'correct', 'condition', 'cond_categoryNontemplate'], 
                header=['onset', 'duration', 'trial_type', 'response_time', 'correct', 'condition', 'condition_nonTTemplate'], index=False)
            # Localizer run (still to do)
            pass

            ###### note that there is an extra column in pp1-6 and the subjects are not ordered!
            ###### THE ONSET times OF THE FIRST SIX SUBJECTS ARE DIFFERENT compared to the rest!!!

            ### Other code (still to do)
            pass