#!/Users/abrahamn/anaconda2/bin/python

import os
from subprocess import call
import numpy as np
import sys

def printbig(s):
        n = len(s)
        print('')
        print('*'*(n+4))
        print('* '+s+' *')
        print('*'*(n+4))
        print('')

run = sys.argv[1]
labels = sys.argv[2:]
print(labels)
analysis_name = '_'.join([l.split('o')[1] for l in labels])
print(analysis_name)

override_swarp = True
override_source_lists = True
override_gen_masks = True
override_coarse_binning = True

raw_frames_root = '/Volumes/abraham/xcor_data/'+run+'/'
analysis_root = '/Volumes/abraham/xcor_data/analysis/'+run+'/' + analysis_name + '/'

ra_cent_deg = 0.
dec_cent_deg = -27

fine_pixel_asec = 1.86
fov_deg = 12
n_fine = int(fov_deg*3600./fine_pixel_asec)

target_coarse_res_asec = 7.*60
coarse_bin_factor = int(np.round(target_coarse_res_asec/fine_pixel_asec))

# create analysis_root if doesn't exist
if not os.path.isdir(analysis_root):
        printbig('making directory'+analysis_root)
        os.makedirs(analysis_root)

# unpack .fits.fz to .fits if not already done
# for l in labels:
# 	if not os.path.isfile(analysis_root + l + '_original.fits'):
# 		printbig('unpacking '+l+'.fits into '+l+'.fits.fz')
# 		print(raw_frames_root+l+'.fits.fz ==> ' + analysis_root+l+'_original.fits')
# 		call(['funpack','-E','1','-O',analysis_root+l+'_original.fits',raw_frames_root+l+'.fits.fz'])

# swarp
if not os.path.isfile(analysis_root + analysis_name + '.fits') or override_swarp:
	call(['swarp']+[raw_frames_root+l+'.fits' for l in labels]+['-c','grid_to_ortho.swarp', \
	'-IMAGEOUT_NAME',analysis_root+analysis_name+'.fits', \
	'-WEIGHTOUT_NAME','/dev/null', \
	'-RESAMPLING_TYPE','NEAREST', '-PIXEL_SCALE',str(fine_pixel_asec), \
	'-IMAGE_SIZE',str(n_fine)+','+str(n_fine),
	'-CENTER',str(ra_cent_deg)+','+str(dec_cent_deg)])


# generate source lists for masking
if not (os.path.isfile(analysis_root+analysis_name+'_artifacts_sources.reg') and \
           os.path.isfile(analysis_root+analysis_name+'_artifacts_sources.csv') and \
           os.path.isfile(analysis_root+analysis_name+'_artifacts.csv')) or override_source_lists:

        printbig('generating source lists for masking')
        call(['./catalog_to_mask_lists.py',analysis_root+analysis_name+'.fits',\
              raw_frames_root+labels[0]+'.dph',analysis_root,labels[0]])

# make fine res fits image of masks for debugging
# and with just bright pixel artifacts masked
if not (os.path.isfile(analysis_root+analysis_name+'_artifacts_sources.fits') and \
          os.path.isfile(analysis_root+analysis_name+'_artifacts.fits')) or override_gen_masks:

        printbig('making fits images of masks')
        call(['./maskcsv2fits.py',analysis_root+analysis_name+'.fits', \
              analysis_root+labels[0]+'_artifacts_sources.csv', \
              analysis_root+analysis_name+'_artifacts_sources.fits'])

        call(['./maskcsv2fits.py',analysis_root+analysis_name+'.fits', \
              analysis_root+labels[0]+'_artifacts.csv', \
              analysis_root+analysis_name+'_artifacts.fits'])

# bin to 3' resolution, excluding masked regions
if not (os.path.isfile(analysis_root+analysis_name+'_mask_artifacts_coarse.fits') and \
        os.path.isfile(analysis_root+analysis_name+'_mask_artifacts_sources_coarse.fits')) or override_coarse_binning:

        printbig('coarse binning by factor of '+str(coarse_bin_factor)+' to'+str(fine_pixel_asec*coarse_bin_factor/60.)+'\' res')
        call(['./coarse_bin_fits_image_with_mask.py',analysis_root+analysis_name+'.fits', \
             analysis_root+analysis_name+'_artifacts_sources.fits', \
             analysis_root+analysis_name+'_mask_artifacts_sources_coarse.fits', \
             analysis_root+analysis_name+'_mask_artifacts_sources_coarse_weights.fits',str(coarse_bin_factor)])

        call(['./coarse_bin_fits_image_with_mask.py',analysis_root+analysis_name+'.fits', \
             analysis_root+analysis_name+'_artifacts.fits', \
             analysis_root+analysis_name+'_mask_artifacts_coarse.fits', \
             analysis_root+analysis_name+'_mask_artifacts_coarse_weights.fits',str(coarse_bin_factor)])

