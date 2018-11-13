
def get_niftis(subject_id, data_dir, task=[], space='mni'):
    """get_niftis gets niftis from a data_dir
    
    Parameters
    ----------
    subject_id : TYPE
        Description
    data_dir : TYPE
        Description
    task : list, optional
        Description
    space : str, optional
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    from bids.grabbids import BIDSLayout # renamed into bids.layout

    layout = BIDSLayout(data_dir)

    nii_files = layout.get(subject=subject_id,
                           task=task,
                           modality="func",
                           extensions=['.nii.gz'])
    # select the proper space
    if space == 'mni':
        nii_files = [
            mf for mf in nii_files if 'space-MNI152NLin2009cAsym_preproc' in mf.filename]
    elif space == 'T1w':
        nii_files = [
            mf for mf in nii_files if 'space-T1w_preproc' in mf.filename]

    nii_files = [nf.filename for nf in nii_files]
    return nii_files


def get_confounds(subject_id, data_dir, task=[]):
    """get_niftis gets niftis from a data_dir
    
    Parameters
    ----------
    subject_id : TYPE
        Description
    data_dir : TYPE
        Description
    task : list, optional
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    from bids.grabbids import BIDSLayout

    layout = BIDSLayout(data_dir)

    confounds_tsv_files = layout.get(subject=subject_id,
                                     task=task,
                                     modality="func",
                                     type="confounds",
                                     extensions=['.tsv'])
    confounds_tsv_files = [cv.filename for cv in confounds_tsv_files]
    return confounds_tsv_files


def get_events(subject_id, data_dir, task=[]):
    """get_niftis gets niftis from a data_dir
    
    Parameters
    ----------
    subject_id : TYPE
        Description
    data_dir : TYPE
        Description
    task : list, optional
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    from bids.grabbids import BIDSLayout

    layout = BIDSLayout(data_dir)

    event_files = layout.get(subject=subject_id,
                             task=task,
                             modality="func",
                             type="events",
                                  extensions=['.tsv'])

    event_files = [ev.filename for ev in event_files]
    return event_files


def get_masks(mask_directory):
    """Summary
    
    Parameters
    ----------
    mask_directory : TYPE
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    import os
    import glob

    mask_files = glob.glob(os.path.join(mask_directory, '*.nii.gz'))
    print(mask_files)
    return mask_files


def resample_rois(mni_roi_files, mni_epi_space_file):
    """Summary
    
    Parameters
    ----------
    mni_roi_files : TYPE
        Description
    mni_epi_space_file : TYPE
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    import os
    import tempfile
    from nilearn.image import resample_to_img

    output_roi_files = []
    for roi in mni_roi_files:
        opf_roi = resample_to_img(roi, mni_epi_space_file)
        opf_roi_fn = os.path.join(tempfile.mkdtemp(), os.path.split(
            roi.replace('.nii.gz', '_epi.nii.gz'))[1])
        opf_roi.to_filename(opf_roi_fn)
        output_roi_files.append(opf_roi_fn)

    return output_roi_files


def create_VWM_anti_pp_workflow(analysis_info, name='VWM-anti'):
    """Summary
    
    Parameters
    ----------
    analysis_info : TYPE
        Description
    name : str, optional
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    import os.path as op
    import nipype.pipeline as pe
    import tempfile
    import glob
    from nipype.interfaces import fsl
    from nipype.interfaces.utility import Function, Merge, IdentityInterface
    from nipype.interfaces.io import SelectFiles, DataSink
    from nipype.interfaces.ants import ApplyTransforms
    from bids.grabbids import BIDSLayout

    # Importing of custom nodes from spynoza packages; assumes that spynoza is installed:
    # pip install git+https://github.com/spinoza-centre/spynoza.git@develop
    from spynoza.filtering.nodes import Savgol_filter, Savgol_filter_confounds
    from spynoza.conversion.nodes import psc
    from spynoza.utils import get_scaninfo, pickfirst
    from utils import mask_nii_2_hdf5, nistats_confound_glm, mask_to_tsv

    input_node = pe.Node(IdentityInterface(
        fields=['bids_directory', 'fmriprep_directory', 'output_directory', 'mask_directory', 'sub_id']), name='inputspec')

    BIDSNiiGrabber = pe.Node(Function(function=get_niftis, input_names=["subject_id",
                                                                        "data_dir", "task", "space"],
                                      output_names=["nii_files"]), name="BIDSNiiGrabber")
    BIDSNiiGrabber.inputs.space = 'mni'

    BIDSEventsGrabber = pe.Node(Function(function=get_events, input_names=["subject_id",
                                                                           "data_dir", "task"],
                                         output_names=["event_files"]), name="BIDSEventsGrabber")
    
    BIDSConfoundsGrabber = pe.Node(Function(function=get_confounds, input_names=["subject_id",
                                                                                 "data_dir", "task"],
                                            output_names=["confounds_tsv_files"]), name="BIDSConfoundsGrabber")
    
    MaskGrabber = pe.Node(Function(function=get_masks, input_names=["mask_directory"],
                                   output_names=["mask_files"]), name="MaskGrabber")

    HDF5PSCMasker = pe.Node(Function(input_names=['in_files', 'mask_files', 'hdf5_file', 'folder_alias'], output_names=['hdf5_file'],
                                     function=mask_nii_2_hdf5),
                            name='hdf5_psc_masker')
    HDF5PSCMasker.inputs.folder_alias = 'psc'
    HDF5PSCMasker.inputs.hdf5_file = op.join(tempfile.mkdtemp(), 'roi.h5')

    HDF5PSCNuisMasker = pe.Node(Function(input_names=['in_files', 'mask_files', 'hdf5_file', 'folder_alias'], output_names=['hdf5_file'],
                                         function=mask_nii_2_hdf5),
                                name='hdf5_psc_nuis_masker')
    HDF5PSCNuisMasker.inputs.folder_alias = 'psc_nuis'

    # HDF5StatsMasker = pe.Node(Function(input_names = ['in_files', 'mask_files', 'hdf5_file', 'folder_alias'], output_names = ['hdf5_file'],
    #                             function = mask_nii_2_hdf5),
    #                             name = 'hdf5_stats_masker')
    # HDF5StatsMasker.inputs.folder_alias = 'stats'

    HDF5ROIMasker = pe.Node(Function(input_names=['in_files', 'mask_files', 'hdf5_file', 'folder_alias'], output_names=['hdf5_file'],
                                     function=mask_nii_2_hdf5),
                            name='hdf5_roi_masker')
    HDF5ROIMasker.inputs.folder_alias = 'rois'

    ConfoundGLM = pe.MapNode(Function(input_names=['nifti_file', 'confounds_file', 'which_confounds'], output_names=['output_pdf', 'output_nifti'],
                                      function=nistats_confound_glm),
                             name='nistats_confound_glm', iterfield=["nifti_file", "confounds_file"])
    ConfoundGLM.inputs.which_confounds = analysis_info['nuisance_columns']

    # VolTransNode = pe.MapNode(interface=fsl.preprocess.ApplyXFM(apply_xfm=False, apply_isoxfm=True, interp='sinc'),
    #                                                     name='vol_trans', iterfield = ['in_file'])

    # VolTransNode = pe.MapNode(interface=ApplyTransforms(transforms='identity', interpolation='LanczosWindowedSinc'),
    #                                                     name='vol_trans', iterfield = ['input_image'])

    ThreshNode = pe.MapNode(fsl.Threshold(thresh=analysis_info['MNI_mask_threshold'], args='-bin', output_datatype='int'),
                            name='thresh', iterfield=['in_file'])

    TSVMasker = pe.MapNode(Function(input_names=['in_file', 'mask_files'], output_names=['out_file'],
                                 function=mask_to_tsv), iterfield=['in_file'],
                        name='tsv_masker')

    ROIResampler = pe.Node(Function(input_names=['mni_roi_files', 'mni_epi_space_file'], output_names=['output_roi_files'],
                                    function=resample_rois),
                           name='roi_resampler')

    sgfilter = pe.MapNode(interface=Savgol_filter,
                          name='sgfilter',
                          iterfield=['in_file'])
    sgfilter_confounds = pe.MapNode(interface=Savgol_filter_confounds,
                                    name='sgfilter_confounds',
                                    iterfield=['confounds'])

    # Both fmri data and nuisances are filtered with identical parameters
    sgfilter.inputs.polyorder = sgfilter_confounds.inputs.polyorder = analysis_info[
        'sgfilter_polyorder']
    sgfilter.inputs.deriv = sgfilter_confounds.inputs.deriv = analysis_info['sgfilter_deriv']
    sgfilter.inputs.window_length = sgfilter_confounds.inputs.window_length = analysis_info[
        'sgfilter_window_length']
    sgfilter.inputs.tr = sgfilter_confounds.inputs.tr = analysis_info['RepetitionTime']

    # set the psc function
    psc.inputs.func = analysis_info['psc_function']

    datasink = pe.Node(DataSink(), name='sinker')
    datasink.inputs.parameterization = False

    ########################################################################################
    # workflow
    ########################################################################################

    # the actual top-level workflow
    VWM_anti_pp_workflow = pe.Workflow(name=name)

    # data source
    VWM_anti_pp_workflow.connect(
        input_node, 'bids_directory', BIDSEventsGrabber, 'data_dir')
    VWM_anti_pp_workflow.connect(input_node, 'sub_id',
                              BIDSEventsGrabber, 'subject_id')
    VWM_anti_pp_workflow.connect(
        input_node, 'fmriprep_directory', BIDSNiiGrabber, 'data_dir')
    VWM_anti_pp_workflow.connect(input_node, 'sub_id',
                              BIDSNiiGrabber, 'subject_id')
    VWM_anti_pp_workflow.connect(
        input_node, 'fmriprep_directory', BIDSConfoundsGrabber, 'data_dir')
    VWM_anti_pp_workflow.connect(input_node, 'sub_id',
                              BIDSConfoundsGrabber, 'subject_id')
    VWM_anti_pp_workflow.connect(
        input_node, 'mask_directory', MaskGrabber, 'mask_directory')

    # filter and psc
    VWM_anti_pp_workflow.connect(BIDSNiiGrabber, 'nii_files', sgfilter, 'in_file')
    VWM_anti_pp_workflow.connect(sgfilter, 'out_file', psc, 'in_file')
    # do the same filtering on confounds
    VWM_anti_pp_workflow.connect(BIDSConfoundsGrabber, 'confounds_tsv_files', sgfilter_confounds, 'confounds')

    # cleanup GLM
    VWM_anti_pp_workflow.connect(psc, 'out_file', ConfoundGLM, 'nifti_file')
    VWM_anti_pp_workflow.connect(
        sgfilter_confounds, 'out_file', ConfoundGLM, 'confounds_file')

    # preparing masks, ANTS and fsl not working correctly
    # ANTs
    # pearl_pp_workflow.connect(BIDSNiiGrabber, ('nii_files', pickfirst), VolTransNode, 'reference_image')
    # pearl_pp_workflow.connect(MaskGrabber, 'mask_files', VolTransNode, 'input_image')
    # fsl
    # pearl_pp_workflow.connect(BIDSNiiGrabber, ('nii_files', pickfirst), VolTransNode, 'reference')
    # pearl_pp_workflow.connect(MaskGrabber, 'mask_files', VolTransNode, 'in_file')
    # pearl_pp_workflow.connect(VolTransNode, 'output_image', ThreshNode, 'in_file')

    VWM_anti_pp_workflow.connect(
        BIDSNiiGrabber, ('nii_files', pickfirst), ROIResampler, 'mni_epi_space_file')
    VWM_anti_pp_workflow.connect(
        MaskGrabber, 'mask_files', ROIResampler, 'mni_roi_files')
    VWM_anti_pp_workflow.connect(
        ROIResampler, 'output_roi_files', ThreshNode, 'in_file')

    # masking data
    VWM_anti_pp_workflow.connect(psc, 'out_file', HDF5PSCMasker, 'in_files')
    VWM_anti_pp_workflow.connect(ThreshNode, 'out_file',
                              HDF5PSCMasker, 'mask_files')

    VWM_anti_pp_workflow.connect(
        ConfoundGLM, 'output_nifti', HDF5PSCNuisMasker, 'in_files')
    VWM_anti_pp_workflow.connect(ThreshNode, 'out_file',
                              HDF5PSCNuisMasker, 'mask_files')
    VWM_anti_pp_workflow.connect(
        HDF5PSCMasker, 'hdf5_file', HDF5PSCNuisMasker, 'hdf5_file')

    # needs stats before we do a masker....
    # pearl_pp_workflow.connect(VolTransNode, 'out_file', HDF5StatsMasker, 'in_files')
    # pearl_pp_workflow.connect(ThreshNode, 'out_file', HDF5StatsMasker, 'mask_files')
    # pearl_pp_workflow.connect(HDF5PSCNuisMasker, 'hdf5_file', HDF5StatsMasker, 'hdf5_file')

    VWM_anti_pp_workflow.connect(
        ROIResampler, 'output_roi_files', HDF5ROIMasker, 'in_files')
    VWM_anti_pp_workflow.connect(ThreshNode, 'out_file',
                              HDF5ROIMasker, 'mask_files')
    VWM_anti_pp_workflow.connect(
        HDF5PSCNuisMasker, 'hdf5_file', HDF5ROIMasker, 'hdf5_file')

    # mask to .tsv, for one timecourse per roi
    VWM_anti_pp_workflow.connect(
        ROIResampler, 'output_roi_files', TSVMasker, 'mask_files')
    VWM_anti_pp_workflow.connect(
        ConfoundGLM, 'output_nifti', TSVMasker, 'in_file')

    # set up output folder
    VWM_anti_pp_workflow.connect(
        input_node, 'output_directory', datasink, 'base_directory')

    # connect all outputs to datasink
    VWM_anti_pp_workflow.connect(
        ConfoundGLM, 'output_nifti', datasink, 'confound_glm')
    VWM_anti_pp_workflow.connect(
        BIDSEventsGrabber, 'event_files', datasink, 'events')
    VWM_anti_pp_workflow.connect(sgfilter, 'out_file', datasink, 'sg_filter')
    VWM_anti_pp_workflow.connect(
        sgfilter_confounds, 'out_file', datasink, 'sg_filter_confound')
    VWM_anti_pp_workflow.connect(psc, 'out_file', datasink, 'psc')
    VWM_anti_pp_workflow.connect(
        ROIResampler, 'output_roi_files', datasink, 'masks_f')
    VWM_anti_pp_workflow.connect(ThreshNode, 'out_file', datasink, 'masks_b')
    VWM_anti_pp_workflow.connect(TSVMasker, 'out_file', datasink, 'tsv')
    VWM_anti_pp_workflow.connect(HDF5PSCNuisMasker, 'hdf5_file', datasink, 'h5')
    VWM_anti_pp_workflow.connect(
        ConfoundGLM, 'output_pdf', datasink, 'confound_glm_report')

    return VWM_anti_pp_workflow
