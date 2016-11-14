#!/Users/abrahamn/anaconda2/bin/python

import commands
from subprocess import call
import numpy as np

for fieldi in range(1):

	labels = open('field'+str(fieldi)+'.labels').read().split('\n')

	# run swarp
	raw_frames_root = '/Volumes/abraham/xcor_data/ATLAS_mwa57694_rereduction/'
	analysis_root = '/Volumes/abraham/xcor_data/analysis/ATLAS_mwa57694_rereduction/field'+str(fieldi)+'/'
	analysis_name = 'field'+str(fieldi)

	ra_cent_deg,dec_cent_deg=np.genfromtxt('field'+str(fieldi)+'.radec',delimiter=',')

	fine_pixel_asec = 1.86
	fov_deg = 5.
	n_fine = int(fov_deg*3600./fine_pixel_asec)

	call(['swarp']+[raw_frames_root+l+'.fits' for l in labels]+\
		['-WEIGHT_TYPE','MAP_WEIGHT']+\
		['-WEIGHT_IMAGE']+[' '.join([raw_frames_root+l+'_weight.fits' for l in labels])]+\
	['-c','stack_at_frame_center.swarp', \
	'-IMAGEOUT_NAME',analysis_root+analysis_name+'.fits', \
	'-WEIGHTOUT_NAME',analysis_root+analysis_name+'_weight.fits', \
	'-RESAMPLING_TYPE','NEAREST', '-PIXEL_SCALE',str(fine_pixel_asec), \
	'-IMAGE_SIZE',str(n_fine)+','+str(n_fine),
	'-CENTER',str(ra_cent_deg)+','+str(dec_cent_deg)])

