def mask_nii_2_hdf5(in_files, mask_files, hdf5_file, folder_alias):
    """masks data in in_files with masks in mask_files,
    to be stored in an hdf5 file

    Takes a list of 3D or 4D fMRI nifti-files and masks the
    data with all masks in the list of nifti-files mask_files.
    These files are assumed to represent the same space, i.e.
    that of the functional acquisitions. 
    These are saved in hdf5_file, in the folder folder_alias.

    Parameters
    ----------
    in_files : list
        list of absolute path to functional nifti-files.
        all nifti files are assumed to have the same ndim
    mask_files : list
        list of absolute path to mask nifti-files.
        mask_files are assumed to be 3D
    hdf5_file : str
        absolute path to hdf5 file.
    folder_alias : str
        name of the to-be-created folder in the hdf5 file.

    Returns
    -------
    hdf5_file : str
        absolute path to hdf5 file.
    """

    import nibabel as nib
    import os.path as op
    import numpy as np
    import tables

    success = True

    mask_data = [np.array(nib.load(mf).get_data(), dtype = bool) for mf in mask_files]
    nifti_data = [nib.load(nf).get_data() for nf in in_files]

    mask_names = [op.split(mf)[-1].split('.nii.gz')[0] for mf in mask_files]
    nifti_names = [op.split(nf)[-1].split('.nii.gz')[0] for nf in in_files]

    h5file = tables.open_file(hdf5_file, mode = "a", title = hdf5_file)
    # get or make group for alias folder
    try:
        folder_alias_run_group = h5file.get_node("/", name = folder_alias, classname='Group')
    except tables.NoSuchNodeError:
        print('Adding group ' + folder_alias + ' to this file')
        folder_alias_run_group = h5file.create_group("/", folder_alias, folder_alias)

    for (roi, roi_name) in zip(mask_data, mask_names):
        # get or make group for alias/roi
        try:
            run_group = h5file.get_node(where = "/" + folder_alias, name = roi_name, classname='Group')
        except tables.NoSuchNodeError:
            print('Adding group ' + folder_alias + '_' + roi_name + ' to this file')
            run_group = h5file.create_group("/" + folder_alias, roi_name, folder_alias + '_' + roi_name)

        h5file.create_array(run_group, roi_name, roi, roi_name + ' mask file for reconstituting nii data from masked data')

        for (nii_d, nii_name) in zip(nifti_data, nifti_names):
            print('roi: %s, nifti: %s'%(roi_name, nii_name))
            n_dims = len(nii_d.shape)
            if n_dims == 3:
                these_roi_data = nii_d[roi]
            elif n_dims == 4:   # timeseries data, last dimension is time.
                these_roi_data = nii_d[roi,:]
            else:
                print("n_dims in data {nifti} do not fit with mask".format(nii_name))
                success = False

            h5file.create_array(run_group, nii_name, these_roi_data, roi_name + ' data from ' + nii_name)

    h5file.close()

    return hdf5_file

def roi_data_from_hdf(data_types_wildcards, roi_name_wildcard, hdf5_file, folder_alias):
    """takes data_type data from masks stored in hdf5_file

    Takes a list of 4D fMRI nifti-files and masks the
    data with all masks in the list of nifti-files mask_files.
    These files are assumed to represent the same space, i.e.
    that of the functional acquisitions. 
    These are saved in hdf5_file, in the folder folder_alias.

    Parameters
    ----------
    data_types_wildcards : list
        list of data types to be loaded.
        correspond to nifti_names in mask_2_hdf5
    roi_name_wildcard : str
        wildcard for masks. 
        corresponds to mask_name in mask_2_hdf5.
    hdf5_file : str
        absolute path to hdf5 file.
    folder_alias : str
        name of the folder in the hdf5 file from which data
        should be loaded.

    Returns
    -------
    output_data : list
        list of numpy arrays corresponding to data_types and roi_name_wildcards
    """
    import tables
    import itertools
    import fnmatch
    import numpy as np
    from IPython import embed as shell

    h5file = tables.open_file(hdf5_file, mode = "r")

    try:
        folder_alias_run_group = h5file.get_node(where = '/', name = folder_alias, classname='Group')
    except tables.NoSuchNodeError:
        # import actual data
        print('No group ' + folder_alias + ' in this file')
        # return None


    all_roi_names = h5file.list_nodes(where = '/' + folder_alias, classname = 'Group')
    roi_names = [rn._v_name for rn in all_roi_names if roi_name_wildcard in rn._v_name]
    if len(roi_names) == 0:
        print('No rois corresponding to ' + roi_name_wildcard + ' in group ' + folder_alias)
        # return None
    
    data_arrays = []
    for roi_name in roi_names:
        try:
            roi_node = h5file.get_node(where = '/' + folder_alias, name = roi_name, classname='Group')
        except tables.NoSuchNodeError:
            print('No data corresponding to ' + roi_name + ' in group ' + folder_alias)
            pass
        all_data_array_names = h5file.list_nodes(where = '/' + folder_alias + '/' + roi_name)
        data_array_names = [adan._v_name for adan in all_data_array_names]
        selected_data_array_names = list(itertools.chain(*[fnmatch.filter(data_array_names, dtwc) for dtwc in data_types_wildcards]))
        
        # if sort_data_types:
        selected_data_array_names = sorted(selected_data_array_names)
        if len(data_array_names) == 0:
            print('No data corresponding to ' + str(selected_data_array_names) + ' in group /' + folder_alias + '/' + roi_name)
            pass
        else:
            print('Taking data corresponding to ' + str(selected_data_array_names) + ' from group /' + folder_alias + '/' + roi_name)
            data_arrays.append([])
            for dan in selected_data_array_names:
                data_arrays[-1].append(eval('roi_node.__getattr__("' + dan + '").read()'))
            print('Taken data corresponding to ' + str(selected_data_array_names) + ' from group /' + folder_alias + '/' + roi_name)
            data_arrays[-1] = np.hstack(data_arrays[-1]) # stack across timepoints or other values per voxel
            if len(data_arrays[-1].shape) == 1:
                data_arrays[-1] = data_arrays[-1][:,np.newaxis]
    all_roi_data_np = np.vstack(data_arrays)    # stack across regions to create a single array of voxels by values (i.e. timepoints)

    h5file.close()

    return all_roi_data_np

def mask_to_tsv(in_file, mask_files):
    import os
    import pandas as pd
    import nibabel as nb
    import numpy as np
    import nilearn.image as image

    mask_names = []
    for mf in mask_files:
        mask_names.append(os.path.split(mf)[1].replace('.nii.gz', ''))

    in_file_nii = nb.load(in_file)
    ipd = in_file_nii.get_data()

    out_file = in_file.replace('.nii.gz', '_roi_ts.tsv')

    # set up output data
    opd = np.zeros((len(mask_files), in_file_nii.shape[-1]))
    # loop across masks and fill out output data
    for mi, maskf in enumerate(mask_files):
        mf_nii = nb.load(maskf)
        md = mf_nii.get_data()
        opd[mi] = np.dot(md.ravel(),ipd.reshape(np.prod(md.shape),-1))/np.sum(md)
    # save out
    opdf = pd.DataFrame(opd.T, columns=mask_names)
    opdf.to_csv(out_file, sep='\t')

    return out_file


def nistats_confound_glm(nifti_file, confounds_file, which_confounds):
    import pandas as pd
    import nibabel as nb
    import numpy as np
    from nistats.regression import OLSModel
    from scipy.stats import zscore
    from nilearn.image import math_img
    from nilearn.plotting import plot_stat_map
    import matplotlib.pyplot as plt

    infile = nb.load(nifti_file)
    mean_img = math_img('np.mean(infile, axis=-1)', infile=infile)
    in_data = infile.get_data().astype(np.float32)

    confounds_table = pd.read_table(confounds_file)[which_confounds]
    # we assume the confounds nor the nifti_file need temporal filtering 
    cfz = confounds_table.apply(zscore)
    cfz['intercept'] = np.ones(infile.shape[-1])

    om = OLSModel(np.array(cfz))
    om_rr = om.fit(in_data.reshape((-1,infile.shape[-1])).T)
    resid_img = nb.Nifti1Image(om_rr.resid.T.reshape(infile.shape).astype(np.float32), affine=infile.affine, header=infile.header)
    cleaned_img = math_img('(resid_img + mean_img[...,np.newaxis]).astype(np.float32)', resid_img=resid_img, mean_img=mean_img)
    output_nifti = nifti_file.replace('.nii.gz', '_nuis.nii.gz')
    cleaned_img.to_filename(output_nifti)
    
    output_pdf = confounds_file.replace('.tsv', '_sd-diff.pdf')
    f = plt.figure(figsize=(24,6))
    plot_stat_map(math_img('(infile.std(axis=-1)-cleaned_img.std(axis=-1))/infile.std(axis=-1)', infile=infile, cleaned_img=cleaned_img), 
        bg_img=mean_img, figure=f, cut_coords=(0,0,0), threshold=0.125, vmax=1, cmap='viridis', output_file=output_pdf)

    return output_pdf, output_nifti

