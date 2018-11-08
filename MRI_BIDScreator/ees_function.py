def _get_pe_index(meta):
	pe = meta['PhaseEncodingDirection']
	try:
		return {'i': 0, 'j': 1, 'k': 2}[pe[0]]
	except KeyError:
		raise RuntimeError('"%s" is an invalid PE string' % pe)


def get_ees(in_meta, in_file=None):
	import nibabel as nb

	# Use case 1: TRT is defined
	trt = in_meta.get('TotalReadoutTime', None)
	if trt is not None:
		return trt

	# All other cases require the parallel acc and npe (N vox in PE dir)
	acc = float(in_meta.get('ParallelReductionFactorInPlane', 1.0))
	npe = nb.load(in_file).shape[_get_pe_index(in_meta)]
	etl = npe // acc

	# Use case 2: TRT is defined
	ees = in_meta.get('EffectiveEchoSpacing', None)
	if ees is not None:
		return ees * (etl - 1)

	# Use case 3 (philips scans)
	wfs = in_meta.get('WaterFatShift', None)
	if wfs is not None:
		fstrength = in_meta['MagneticFieldStrength']
		wfd_ppm = 3.4  # water-fat diff in ppm
		g_ratio_mhz_t = 42.57  # gyromagnetic ratio for proton (1H) in MHz/T
		wfs_hz = fstrength * wfd_ppm * g_ratio_mhz_t
		return wfs / wfs_hz

	raise ValueError('Unknown total-readout time specification')

def readWfs(pathIN):
	import re
	f = open(pathIN, "r")
	contents = f.read()
	idx = contents.index("pixels")
	waterFatShift = [float(s) for s in re.findall(r'-?\d+\.?\d*',contents[idx+22:idx+32])][0]

	return waterFatShift

def calculateParameters(pathIN, waterFatShift = False, effectiveEchoSpacing = False):
	if waterFatShift and effectiveEchoSpacing:
		raise Exception('Do not provide both, Water Fat Shift and Effective Echo Spacing!')
	elif waterFatShift:
		in_meta = {'WaterFatShift': waterFatShift,
				'MagneticFieldStrength': 3,
				'PhaseEncodingDirection': 'j-',
				'ParallelReductionFactorInPlane': 1.5}
	elif effectiveEchoSpacing:
		in_meta = {'EffectiveEchoSpacing': effectiveEchoSpacing,
				'MagneticFieldStrength': 3,
				'PhaseEncodingDirection': 'j-',
				'ParallelReductionFactorInPlane': 1.5}
	else:
		print('Did not calculate anything...')
	return get_ees(in_meta, pathIN)

# For manuell call, not function call
# in_meta = {'WaterFatShift': 13.121,
# 		   'MagneticFieldStrength': 3,
# 		   'PhaseEncodingDirection': 'j-',
# 		   'ParallelReductionFactorInPlane': 1.5}

# # For philips
# in_file = 'sub-01_dir-NR_run-1_epi.nii.gz'
# fmap_file = 'C:\PycharmProjects\FMRI_NRosT\sub-01\Fmap'+'\\'+in_file
# print(get_ees(in_meta, fmap_file))
