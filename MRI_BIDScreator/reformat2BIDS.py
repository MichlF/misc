'''
This script ...
- ...reads in nii files from some folder and reformats them into the BIDS format.
- ...creates JSON files for .epi files.
- ...reads in behavioral .csv files from OpenSesame and creates .tsv event files
        which can be used as design matrices for nistats/nilearn

In the future it will also be able to ...
- ...create dataset_description.json and task jsons.
- ...read in and format eye tracking data.
- ...read in and format physiology data.

@Michel Failing, 2019
'''

# Import
import os, errno, subprocess, re, json, B0_separation, ees_function, shutil
import nibabel as nb
import pandas as pd
from shutil import copyfile
from pprint import pprint

# Some helper functions
def folderCreator(path):
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

### Start script
# FOR NOW:
use_PP = [13] # Only analyse subjects in this list of subjects
sesIdx = '01' # Session (currently only one is working)
nrRuns = 9 # Number of separate runs (8 + 1 localizer here)
no_trialsBlock = 24 # trials per experimental block
time_memEvent = 7.7 # memory event interval (in sec)
time_beginMemEvent = 10.5 # begin of first memory event (in sec)
time_trial = 19.6 # duration of an experimental trial (in sec)
time_trialLoc = .250 # duration of a localizer "trial" (in sec)
time_TR = .7 # duration of one repetition time

# What to do?
parrec2nii = False # Convert .PAR/.REC to nifti files
construct = True # Construct the BIDS structure (includes the moving and renaming of nifti files)
fixNiftiHeader = True # If BIDS validator throws errors because of the TR being defined in ms instead of s
behavior = True # Write .tsv files

# Paths
pathBEHdata = '/Users/michlf/Dropbox/Work/Data/fMRI/NegativeTemplate/forScanner/mri/beh'
pathMRIdata = '/Volumes/VUHDD/NROST_data/' 
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
    # Remove all irrelevant folders
    folderlist = [s for s in folderlist if 'NROST_VWM' in s[:9]]
    folderCreator(pathBIDS+'code/')
    folderCreator(pathBIDS+'derivatives/')

    # Create folder structure for any subject found
    subList = []
    for i in range(len(folderlist)):
        try:
            subList.append(folderlist[i].split("sub-",1)[1])
            nrSubj = folderlist[i].split("sub-",1)[1]
            if use_PP and int(nrSubj) not in use_PP:
                print('\n...skipped subject {0}...\n'.format(nrSubj))
                continue # if not a subject to use when specified, move on to next iteration
            print('\n...running subject {0}...\n'.format(nrSubj))
        except Exception as e:
            print(e)
            continue # if not a proper subject folder, move on to next iteration
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
                        # Although the validator triggers an error, it does not need to be fixed.
                        # Let's do it nevertheless...
                        if fixNiftiHeader: # BIDS validator throws an error because the time units of the repetition time are in ms rather than in s (seconds is standard for the JSON)
                            fixed = nb.load(path2New)
                            fixed.header.set_xyzt_units(8) # It is on 18 for some reason...
                            nb.save(fixed, path2New)
                        # JSON
                        path_fmap_PAR = pathMRIdata + folderlist[i] + '/' + [s for s in os.listdir(pathMRIdata+folderlist[i]) if "{0}_bold_".format(runIdx) in s and ".PAR" in s][0]
                        path_intended = '"ses-{2}/func/sub-{0}_ses-{2}_task-NRoST_run-{1}_bold.nii.gz"'.format(nrSubj, runIdx, sesIdx)
                        path_fmap = path_fmapDir + 'sub-{0}_ses-{2}_dir-NR_run-{1}'.format(nrSubj, runIdx, sesIdx)
                    except Exception as e:
                        print('Error message: ', e)
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
                            fixed.header.set_xyzt_units(8) # It is on 18 for some reason...
                            nb.save(fixed, path2New)
                        # JSON
                        path_fmap_PAR = pathMRIdata + folderlist[i] + '/' + [s for s in os.listdir(pathMRIdata+folderlist[i]) if "{0}-bold_".format('Local') in s and ".PAR" in s][0]
                        path_intended = '"ses-{2}/func/sub-{0}_ses-{2}_task-Localizer_run-{1}_bold.nii.gz"'.format(nrSubj, runIdx, sesIdx)
                        path_fmap = path_fmapDir + 'sub-{0}_ses-{2}_dir-NR_run-{1}'.format(nrSubj, runIdx, sesIdx)
                    except Exception as e:
                        print('Error message: ', e)
                        print('Trying to access different session for localizer...')
                        try:
                            # Create folder for new session
                            folderCreator(pathSubj+'/ses-{0}/anat/'.format('02'))
                            folderCreator(pathSubj+'/ses-{0}/func/'.format('02'))
                            folderCreator(pathSubj+'/ses-{0}/fmap/'.format('02'))
                            folderCreator(pathSubj+'/ses-{0}/beh/'.format('02'))
                            # EPI
                            path2Original = pathMRIdata + folderlist[i] + '/localizer_rerun/' + [s for s in os.listdir(pathMRIdata+folderlist[i] + '/localizer_rerun/') if "{0}_epi_".format('Local') in s and ".nii.gz" in s][0]
                            path2New = pathSubj + '/ses-{2}/fmap/sub-{0}_ses-{2}_dir-NR_run-{1}_epi.nii.gz'.format(nrSubj, runIdx, '02')
                            copyfile(path2Original, path2New)
                            # BOLD
                            path2Original = pathMRIdata + folderlist[i] + '/localizer_rerun/' + [s for s in os.listdir(pathMRIdata+folderlist[i] + '/localizer_rerun/') if "{0}-bold_".format('Local') in s and ".nii.gz" in s][0]
                            path2New = pathSubj + '/ses-{2}/func/sub-{0}_ses-{2}_task-Localizer_run-{1}_bold.nii.gz'.format(nrSubj, runIdx, '02')
                            copyfile(path2Original, path2New)
                            # Fix nifti header
                            if fixNiftiHeader: # BIDS validator throws an error because the time units of the repetition time are in ms rather than in s (seconds is standard for the JSON)
                                fixed = nb.load(path2New)
                                fixed.header.set_xyzt_units(8) # It is on 18 for some reason...
                                nb.save(fixed, path2New)
                            # JSON
                            path_fmap_PAR = pathMRIdata + folderlist[i] + '/localizer_rerun/' + [s for s in os.listdir(pathMRIdata+folderlist[i] + '/localizer_rerun/') if "{0}-bold_".format('Local') in s and ".PAR" in s][0]
                            path_intended = '"ses-{2}/func/sub-{0}_ses-{2}_task-Localizer_run-{1}_bold.nii.gz"'.format(nrSubj, runIdx, '02')
                            path_fmap = pathSubj + '/ses-{2}/fmap/sub-{0}_ses-{2}_dir-NR_run-{1}'.format(nrSubj, runIdx, '02')
                            # B0 separation
                            # If we come till here, there was a localizer rerun. So we will also have to separate B0 for this session as well
                            path2Original = pathMRIdata + folderlist[i] + '/localizer_rerun/' + [s for s in os.listdir(pathMRIdata+folderlist[i] + '/localizer_rerun/') if "B0" in s and ".nii.gz" in s][0]
                            B0_separation.separation(path2Original, pathSubj+'/ses-{0}/fmap/'.format('02'), nrSubj, '02')
                            # Report
                            print('Found localizer. Created second session and did B0 separation.')
                        
                        except Exception as e:
                            print("Error message:  ", e)
                            # Remove unnecessarily created folders
                            shutil.rmtree(pathSubj+'/ses-{0}/anat/'.format('02'))
                            shutil.rmtree(pathSubj+'/ses-{0}/func/'.format('02'))
                            shutil.rmtree(pathSubj+'/ses-{0}/fmap/'.format('02'))
                            shutil.rmtree(pathSubj+'/ses-{0}/beh/'.format('02'))
                            #if os.path.isdir(pathMRIdata + folderlist[i] + '/' + [s for s in os.listdir(pathMRIdata+folderlist[i]) if "{0}_epi_".format('Local') in s and ".nii.gz" in s][0]
                            print('Creating localizer files failed. File missing?')

                ### JSON files ###
                # Individual EPIs (for each experimental run incl. localizer)
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
                    print('Error message: ', e)
                    print('Creating JSON EPIs failed. File missing?')

            # Phasediff (only needed once, so out of the loop)
            File = open(path_fmapDir + 'sub-{0}_ses-{1}_phasediff.json'.format(nrSubj, sesIdx), 'w')
            File.write(''.join('{')+ '\n')
            File.write(''.join('    "EchoTime1": 0.003,')+ '\n')
            File.write(''.join('    "EchoTime2": 0.008,')+ '\n')
            iterRun = nrRuns
            if os.path.isdir(pathSubj+'/ses-{0}'.format('02')):
                iterRun = nrRuns-1
            for run in range(iterRun):
                if run == 0:
                    File.write(''.join('    "IntendedFor": [' + path_intendedList[run] + ',')+ '\n')
                elif run == iterRun-1:
                    File.write(''.join('                    ' + path_intendedList[run] + ']')+ '\n')
                else:
                    File.write(''.join('                    ' + path_intendedList[run] + ',')+ '\n')
            File.write(''.join('}')+ '\n')
            File.close()

            # If we had the case of a localizer rerun, we need to create another phase difference JSON
            if os.path.isdir(pathSubj+'/ses-{0}'.format('02')):
                File = open(pathSubj + '/ses-{1}/fmap/sub-{0}_ses-{1}_phasediff.json'.format(nrSubj, '02'), 'w')
                File.write(''.join('{')+ '\n')
                File.write(''.join('    "EchoTime1": 0.003,')+ '\n')
                File.write(''.join('    "EchoTime2": 0.008,')+ '\n')
                File.write(''.join('    "IntendedFor": [' + path_intendedList[-1] + ']')+ '\n')
                File.write(''.join('}')+ '\n')
                File.close()

        ### Physiology ###
        pass

        ### Eye-tracking ###
        pass

        ### Behavior ###
        if behavior:
            # Experimental run(s)
            path2Original = pathBEHdata + '/' + [s for s in os.listdir(pathBEHdata) if "subject-{0}_".format(int(nrSubj)) in s and ".csv" in s][0]
            df = pd.read_csv(path2Original)
            # Create/change a few variables to comply with BIDS
            df['duration'] = time_memEvent
            df['condition'] = 'memory'
            df['responseTime'] = round(df['responseTime'])
            df['trial_type'] = df['cond_template'].map(str) + '-' + df['cond_category']
            #df['modulation'] = 1

            # FIR ANALYSES
            # For the .tsv file for the FIR analyses we need to expand the dataframe and thus rather build a new one from scratch
            df1 = pd.read_csv(path2Original, usecols=['responseTime','correct','cond_categoryNontemplate', 'cond_template', 'cond_category'])
            df1['trial_type']= df1['cond_template'].map(str) + '-' + df1['cond_category']
            df_temp = pd.DataFrame(columns=list(df1)) # df with only headers
            for i in range(len(df1.index)):
                for k in range(28): # This method is insanely slow.
                    df_temp = df_temp.append(df1.iloc[i], ignore_index = True) # fill rows with TR times the same info
                    #FIR
                    df_temp.iloc[i*28+k, df_temp.columns.get_loc('trial_type')] = df1.iloc[i, df1.columns.get_loc('trial_type')] + '-{:02d}'.format(k+1) # trial_type label should reflect TR for a given condition
                    #FIR2
                    df_temp.iloc[i*28+k, df_temp.columns.get_loc('trial_type')] = 'TR-{:02d}'.format(k+1) # trial_type label should reflect TR
            # Fill up the rest
            df_temp['duration'] = time_TR
            df_temp['responseTime'] = round(df['responseTime'])
            #df_temp['modulation'] = 1
            
            # Now, chop the file up into individual runs
            for runIdx in range(1,nrRuns):
                # Cut the runs
                df_final = df.truncate(before=(runIdx-1)*no_trialsBlock, after=no_trialsBlock-1+((runIdx-1)*no_trialsBlock))
                df_final_FIR = df_temp.truncate(before=(runIdx-1)*no_trialsBlock*28, after=no_trialsBlock*28-1+((runIdx-1)*no_trialsBlock*28))
                if int(nrSubj) < 7: # One TR at the end of the first trial of each run is missing
                    df_final['onset'] = [round(num * time_trial + time_beginMemEvent - time_TR, 2) for num in range(no_trialsBlock)] # adjust all onsets by removing 1 TR
                    df_final.iloc[0, df_final.columns.get_loc('onset')] = 10.5 # change the first onset
                    # FIR model: We want each of the 28 TRs per trial to be a regressor
                    df_final_FIR['onset'] = [round(num * time_TR + time_beginMemEvent, 2) for num in range(no_trialsBlock*28)] # adjust all onsets by removing 1 TR
                    df_final_FIR.iloc[0, df_final_FIR.columns.get_loc('onset')] = 10.5 # change the first onset
                else:
                    df_final['onset'] = [round(num * time_trial + time_beginMemEvent, 2) for num in range(no_trialsBlock)]
                    df_final_FIR['onset'] = [round(num * time_TR + time_beginMemEvent, 2) for num in range(no_trialsBlock*28)]
                # Save
                # Normal
                path2New = pathSubj + '/ses-{2}/func/sub-{0}_ses-{2}_task-NRoST_run-{1}_events.tsv'.format(nrSubj, runIdx, sesIdx)
                df_final.to_csv(path2New, sep='\t', columns=['onset', 'duration', 'trial_type', 'responseTime', 'correct', 'condition', 'cond_categoryNontemplate'], 
                header=['onset', 'duration', 'trial_type', 'response_time', 'correct', 'condition', 'condition_nonTTemplate'], index=False)
                # FIR
                path2New = pathSubj + '/ses-{2}/func/sub-{0}_ses-{2}_task-NRoST_run-{1}_events_FIR2.tsv'.format(nrSubj, runIdx, sesIdx)
                df_final_FIR.to_csv(path2New, sep='\t', columns=['onset', 'duration', 'trial_type', 'responseTime', 'correct', 'cond_template', 'cond_category', 'cond_categoryNontemplate'], 
                header=['onset', 'duration', 'trial_type', 'response_time', 'correct', 'cond_template', 'cond_category', 'condition_nonTTemplate'], index=False)

            # Localizer(s)
            try:
                path2Original = pathBEHdata + '/' + [s for s in os.listdir(pathBEHdata) if "subject-{0}loc_".format(int(nrSubj)) in s and ".csv" in s][0]
            except:
                print('WARNING: could not find .csv for localizer. Skipped this subject!')
                continue
            df = pd.read_csv(path2Original)
            # Construct a dataframe for the events.tsv file
            df_final = pd.DataFrame(columns=['onset','duration','trial_type','exemplarOrderFulllist', 'cond_order', 'block_list'])
            if (int(nrSubj)<17 and int(nrSubj)!=14) or int(nrSubj)==18: # convenience function logged incorrectly in these subjects, so we will have to do it manually
                # Clean variable strings
                block_list = re.sub("[']", '', df['block_list'][0])
                block_list = block_list[1:-1].split(", ")
                while 'fix' in block_list: block_list.remove('fix')
                cond_order = df['cond_order'][0].split(", ")
                exemplar_order = re.sub("[']", '', df['exemplarOrderFulllist'][0])
                exemplar_order = exemplar_order[1:-1].split(", ")
                while '99' in exemplar_order: exemplar_order.remove('99')
                while '[99' in exemplar_order: exemplar_order.remove('[99')
                while '99]' in exemplar_order: exemplar_order.remove('99]')

                # Loop through the trial to create master string and then construct the vector
                idx,idx2,trialTypes = 0,0,[]
                for block in range(len(block_list)):
                    strCond = block_list[block]
                    for cond in range(2):
                        if 'intact' in cond_order[idx]:
                            strBlock = 'Real'
                        else:
                            strBlock = 'Scram'
                        idx+=1
                        for trial in range(20):
                            exemplar_order[idx2] = exemplar_order[idx2].replace('[','')
                            exemplar_order[idx2] = exemplar_order[idx2].replace(']','')
                            strTrial = str(int(exemplar_order[idx2])+1)
                            trialTypes.append(strCond+strBlock+strTrial)
                            idx2+=1
                # Now add fixation sequences
                trialTypes = 40*['fix']+trialTypes[:120]+40*['fix']+trialTypes[120:240]+40*['fix']+trialTypes[240:360]+40*['fix']+trialTypes[360:480]
            else:
                trialTypes = re.sub("[']", '', df['pathsFulllist'][0])
                trialTypes = trialTypes[1:-1].split(", ")
            df_final['onset'] = [num * time_trialLoc*2 + time_beginMemEvent for num in range(len(trialTypes))]
            df_final['duration'] = time_trialLoc
            df_final['trial_type'] = [trialTypes[trial] for trial in range(len(trialTypes))]
            #df_final['modulation'] = 1
            # Save
            if os.path.isdir(pathSubj+'/ses-{0}'.format('02')):
                path2New = pathSubj + '/ses-{2}/func/sub-{0}_ses-{2}_task-Localizer_run-{1}_events.tsv'.format(nrSubj, nrRuns, '02')
            else:
                path2New = pathSubj + '/ses-{2}/func/sub-{0}_ses-{2}_task-Localizer_run-{1}_events.tsv'.format(nrSubj, nrRuns, sesIdx)
            df_final.to_csv(path2New, sep='\t', columns=['onset', 'duration', 'trial_type'], index=False)

#END
