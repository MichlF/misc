'''
This script ...
- ...reads in nii files from some folder and reformats them into the BIDS format.
- ...creates JSON files for .epi files.
- ...reads in behavioral .csv files from OpenSesame and creates .tsv event files.

In the future it will also be able to ...
- ...create dataset_description.json and task jsons.
- ...read in and format eye tracking data.
- ...read in and format physiology data.

@Michel Failing, 2018
'''

# Import
import os, subprocess, re, json, B0_separation, ees_function
import nibabel as nb
import pandas as pd
from shutil import copyfile
from pprint import pprint

# Some helper functions
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
time_trialLoc = .250

# What to do?
parrec2nii = False
construct = False
fixNiftiHeader = False
behavior = True

# Paths
pathBEHdata = '/Users/michlf/Dropbox/Work/Data/fMRI/NegativeTemplate/forScanner/mri/beh'
pathMRIdata = '/Users/michlf/NROST_working/' 
pathBIDS = '/Users/michlf/NROST_working/fMRI_NROST/' #'/Users/michlf/Documents/GitHub/fMRI_NRoST

# Start reformatting
# Implement parrec2nii (still to do...)
if parrec2nii:
    pass
    #p = subprocess.Popen(['parrec2nii', '-c', '*.PAR'], cwd=pathMRIdata)
    #os.system("cd {0}".format(pathMRIdata))
    #os.system('ls')
    #os.system('parrec2nii -c *.PAR')

# Reformat
if any('sub-' in s for s in os.listdir(pathMRIdata)):  # check whether there is any subject data

    # Some other preparations
    os.chdir(pathMRIdata)
    folderlist = os.listdir(pathMRIdata)
    folderlist.sort()
    folderCreator(pathBIDS+'code/')
    folderCreator(pathBIDS+'derivatives/')

    # Create folder structure for any subject found
    subList = []
    for i in range(len(folderlist)):
        try:
            subList.append(folderlist[i].split("sub-",1)[1])
            nrSubj = folderlist[i].split("sub-",1)[1]
            print('\n...running subject {0}...\n'.format(nrSubj))
        except:
            continue # if not a subject folder, move on to next iteration
        pathSubj = pathBIDS + 'sub-' + nrSubj
        folderCreator(pathSubj+'/ses-{0}/anat/'.format(sesIdx))
        folderCreator(pathSubj+'/ses-{0}/func/'.format(sesIdx))
        folderCreator(pathSubj+'/ses-{0}/fmap/'.format(sesIdx))
        folderCreator(pathSubj+'/ses-{0}/beh/'.format(sesIdx))

        ### Construct name, then copy and rename files ###
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
            B0_separation.separation(path2Original, pathSubj+'/ses-{0}/fmap/'.format(sesIdx), nrSubj, sesIdx)
            
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
                        # Fix nifti header
                        if fixNiftiHeader: # BIDS validator throws an error because the time units of the repetition time are in ms rather than in s (seconds is standard for the JSON)
                            fixed = nb.load(path2New)
                            fixed.header.set_xyzt_units(8) # It is on 18...
                            nb.save(fixed, path2New)
                        # JSON
                        path_fmap_PAR = pathMRIdata + folderlist[i] + '/' + [s for s in os.listdir(pathMRIdata+folderlist[i]) if "{0}_bold_".format(runIdx) in s and ".PAR" in s][0]
                        path_intended = '"ses-{2}/func/sub-{0}_ses-{2}_task-NRoST_run-{1}_bold.nii.gz"'.format(nrSubj, runIdx, sesIdx)
                        path_fmap = path_fmapDir + 'sub-{0}_ses-{2}_dir-NR_run-{1}'.format(nrSubj, runIdx, sesIdx)
                    except Exception as e:
                        print('Error message:  ', e)
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
                        # Fix nifti header
                        if fixNiftiHeader: # BIDS validator throws an error because the time units of the repetition time are in ms rather than in s (seconds is standard for the JSON)
                            fixed = nb.load(path2New)
                            fixed.header.set_xyzt_units(8) # It is on 18...
                            nb.save(fixed, path2New)
                        # JSON
                        path_fmap_PAR = pathMRIdata + folderlist[i] + '/' + [s for s in os.listdir(pathMRIdata+folderlist[i]) if "{0}-bold_".format('Local') in s and ".PAR" in s][0]
                        path_intended = '"ses-{2}/func/sub-{0}_ses-{2}_task-Localizer_run-{1}_bold.nii.gz"'.format(nrSubj, runIdx, sesIdx)
                        path_fmap = path_fmapDir + 'sub-{0}_ses-{2}_dir-NR_run-{1}'.format(nrSubj, runIdx, sesIdx)
                    except Exception as e:
                        print('Error message:  ', e)
                        print('Localizer run did not work. File missing?')

            ### JSON files ###
                # Individual EPIs (for each run incl. localizer)
                try:
                    wfs = ees_function.readWfs(path_fmap_PAR)
                    print('wfs for run {0}: {1} px'.format(runIdx, wfs))
                    EES = ees_function.calculateParameters(path_fmap+'_epi.nii.gz', waterFatShift=wfs)
                    totalReadoutTime = ees_function.calculateParameters(path_fmap+'_epi.nii.gz', effectiveEchoSpacing=EES)
                    File = open(path_fmap + '_epi.json', 'w')
                    File.write(''.join('{')+ '\n')
                    File.write(''.join('    "PhaseEncodingDirection": "j-",')+ '\n')
                    File.write(''.join('    "EffectiveEchoSpacing": ' + '{0}'.format(EES) +',')+ '\n')
                    File.write(''.join('    "TotalReadoutTime": ' + '{0}'.format(totalReadoutTime) +',')+ '\n')
                    File.write(''.join('    "IntendedFor": [' + path_intended + ']')+ '\n')
                    File.write(''.join('}')+ '\n')
                    File.close()
                    path_intendedList.append(path_intended) # save for phase difference JSON
                except Exception as e:
                    print('Error message:  ', e)
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

        ### Physiology (still to do) ###
        pass

        ### Behavior ###
        if behavior:
            # Experimental run(s)
            path2Original = pathBEHdata + '/' + [s for s in os.listdir(pathBEHdata) if "subject-{0}_".format(int(nrSubj)) in s and ".csv" in s][0]
            df = pd.read_csv(path2Original)
            # Create/change a few variables to comply with BIDS
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
                    df_final.iloc[0, df_final.columns.get_loc('onset')] = 10.5 # change the first onset
                else:
                    df_final['onset'] = [round(num * time_trial + time_beginMemEvent, 2) for num in range(no_trialsBlock)]
                # Save
                path2New = pathSubj + '/ses-{2}/func/sub-{0}_ses-{2}_task-NRoST_run-{1}_events.tsv'.format(nrSubj, runIdx, sesIdx)
                df_final.to_csv(path2New, sep='\t', columns=['onset', 'duration', 'trial_type', 'responseTime', 'correct', 'condition', 'cond_categoryNontemplate'], 
                header=['onset', 'duration', 'trial_type', 'response_time', 'correct', 'condition', 'condition_nonTTemplate'], index=False)

            # Localizer(s)
            path2Original = pathBEHdata + '/' + [s for s in os.listdir(pathBEHdata) if "subject-{0}loc_".format(int(nrSubj)) in s and ".csv" in s][0]
            df = pd.read_csv(path2Original)
            # Construct a dataframe for the events.tsv file
            df_final = pd.DataFrame(columns=['onset','duration','trial_type'])
            trialTypes = re.sub("[']", '', df['pathsFulllist'][0])
            trialTypes = trialTypes[1:-1].split(", ")
            test=re.sub("[']", '', df['exemplarOrderFulllist'][0])
            test=test[1:-1].split(", ")
            print(test, len(test), type(test), test[40], 'bla', test[39])
            df_final['onset'] = [num * time_trialLoc*2 + time_beginMemEvent for num in range(len(trialTypes))]
            df_final['duration'] = time_trialLoc
            df_final['trial_type'] = [trialTypes[trial] for trial in range(len(trialTypes))]
            # Save
            path2New = pathSubj + '/ses-{2}/func/sub-{0}_ses-{2}_task-Localizer_run-{1}_events.tsv'.format(nrSubj, nrRuns, sesIdx)
            df_final.to_csv(path2New, sep='\t', columns=['onset', 'duration', 'trial_type'], index=False)
