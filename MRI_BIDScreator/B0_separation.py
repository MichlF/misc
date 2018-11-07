'''
This function separates the B0 scan into magnitude and phase images
'''

def separation(pathIN, pathOUT, nrSubj=999, sessID=''):
    import nibabel as nb

    # Separate
    bf = nb.load(pathIN)
    h = bf.header
    af = bf.affine
    f1 = nb.Nifti1Image(bf.get_data()[..., 0], header=h, affine=af)
    f1.to_filename(pathOUT + 'sub-{0}_ses-{1}_magnitude1.nii.gz'.format(nrSubj,sessID))
    f2 = nb.Nifti1Image(bf.get_data()[..., 1], header=h, affine=af)
    f2.to_filename(pathOUT + 'sub-{0}_ses-{1}_phasediff.nii.gz'.format(nrSubj,sessID))
