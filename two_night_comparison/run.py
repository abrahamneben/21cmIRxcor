#!/Users/abrahamn/anaconda2/bin/python

from subprocess import call
import numpy as np

labels = [open('night1.labels').read().split('\n'),open('night2.labels').read().split('\n')]

fine_pixel_asec = 1.86
fov_deg = 3.
n_fine = int(fov_deg*3600./fine_pixel_asec)

raw_frames_roots = ['/Volumes/abraham/xcor_data/ATLAS_mwa57639/','/Volumes/abraham/xcor_data/ATLAS_mwa57694_rereduction/']
analysis_root = '/Volumes/abraham/xcor_data/analysis/two_night_comparison/'

ra_cent_deg = -3.4
dec_cent_deg = -23.5

target_coarse_res_asec = 6.*60
coarse_bin_factor = int(np.round(target_coarse_res_asec/fine_pixel_asec))

for i in range(2):
	analysis_name = 'night'+str(i+1)

	# stack

	call(['swarp']+[raw_frames_roots[i]+l+'.fits' for l in labels[i]]+\
		['-WEIGHT_TYPE','MAP_WEIGHT']+\
		['-WEIGHT_IMAGE']+[' '.join([raw_frames_roots[i]+l+'_weight.fits' for l in labels[i]])]+\
	['-c','two_night_comparison.swarp', \
	'-IMAGEOUT_NAME',analysis_root+analysis_name+'.fits', \
	'-WEIGHTOUT_NAME',analysis_root+analysis_name+'_weight.fits', \
	'-RESAMPLING_TYPE','NEAREST', '-PIXEL_SCALE',str(fine_pixel_asec), \
	'-IMAGE_SIZE',str(n_fine)+','+str(n_fine),
	'-CENTER',str(ra_cent_deg)+','+str(dec_cent_deg)])

	# generate source lists for masking
	call(['../stack_orthogrid_mask_coarsegrid/catalog_to_mask_lists.py',analysis_root+analysis_name+'.fits',\
             raw_frames_roots[i]+labels[i][0]+'.dph',analysis_root,analysis_name])

	# generate mask images from mask lists
	call(['../stack_orthogrid_mask_coarsegrid/maskcsv2fits.py',analysis_root+analysis_name+'.fits', \
	      analysis_root+analysis_name+'_artifacts_sources.csv', \
	      analysis_root+analysis_name+'_artifacts_sources.fits'])

	# coarse grid the unmasked pixels
	call(['../stack_orthogrid_mask_coarsegrid/coarse_bin_fits_image_with_mask.py',analysis_root+analysis_name+'.fits', \
             analysis_root+analysis_name+'_artifacts_sources.fits', \
             analysis_root+analysis_name+'_mask_artifacts_sources_coarse.fits', \
             analysis_root+analysis_name+'_mask_artifacts_sources_coarse_weights.fits',str(coarse_bin_factor)]) 
