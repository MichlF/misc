# from nipype import config
# config.enable_debug_mode()
# Importing necessary packages
import os
import os.path as op
import json
from nipype import config, logging

from post_fmriprep_workflow import create_VWM_anti_pp_workflow

# export SUBJECTS_DIR=/home/shared/2017/reward/pearl_3T/FS_SJID
# for s in {01..40}
# do
#     echo sub-$s
#     python preprocessing.py sub-$s;
# done


# the subject id and experiment vars are commandline arguments to this script.
for sub_id in ['1']:
    analysis_info = {}
    with open('analysis_parameters.json', 'r') as f:
        json_s = f.read()
        analysis_info.update(json.loads(json_s))

    with open('preprocessing_parameters.json', 'r') as f:
        json_s = f.read()
        analysis_info.update(json.loads(json_s))

    opd = op.join(analysis_info['output_directory'], 'sub-%s'%sub_id)

    # we set up the folders and logging there.
    if not op.isdir(opd):
        try:
            os.makedirs(op.join(opd, 'log'))
        except OSError:
            pass

    config.update_config({  'logging': {
                                        'log_directory': op.join(opd, 'log'),
                                        'log_to_file': True,
                                        'workflow_level': 'DEBUG',
                                        'interface_level': 'DEBUG'
                                      },
                            'docker run -v $HOME:$HOME -p ${external_portnr}:${internal_portnr} --expose=${external_portnr} -v ${data_directory_host}:${data_directory_container} \
-v ${code_directory_host}:${code_directory_container} -v /etc/group:/etc/group:ro -v /etc/passwd:/etc/passwd:ro -u=$UID:$(id -g $USER) -i -t knapenlab/nd:${version}execution': {
                                        'stop_on_first_crash': True
                                        }
                        })
    logging.update_logging(config)

    # the actual workflow
    VWM_anti_pp_workflow = create_VWM_anti_pp_workflow(analysis_info, name = 'VWM-anti')

    # set subject
    VWM_anti_pp_workflow.inputs.inputspec.sub_id = sub_id

    # standard output variables
    VWM_anti_pp_workflow.inputs.inputspec.fmriprep_directory = analysis_info['fmriprep_directory']
    VWM_anti_pp_workflow.inputs.inputspec.bids_directory = analysis_info['bids_directory']
    VWM_anti_pp_workflow.inputs.inputspec.output_directory = opd
    VWM_anti_pp_workflow.inputs.inputspec.mask_directory = analysis_info['mask_directory']

    # write out the graph and run
    # VWM_anti_pp_workflow.write_graph(op.join(opd, 'graph.png'))
    VWM_anti_pp_workflow.run('MultiProc', plugin_args={'n_procs': 8})