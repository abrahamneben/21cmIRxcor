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

overwrite_swarp = False
overwrite_mask_source_lists = True
show_fits_image_with_masked_regions = False
overwrite_mask_fits_images = True
override_coarse_binning = True

if len(sys.argv) != 2:
	print('Usage: ./run.py [02a57639o0338I]')
	sys.exit(0)

#label = '02a57639o0338I'
#label = '02a57639o0347I'
label = sys.argv[1]
raw_frames_root = '/Volumes/abraham/xcor_data/ATLAS_mwa57639/'
analysis_root = '/Volumes/abraham/xcor_data/analysis/' + label + '/'
printbig('processing '+label)

ra_cent_deg = 0
dec_cent_deg = -27

fine_pixel_asec = 1.86 # same as original stereographic projection res
fov_deg = 15.
n_fine = int(fov_deg*3600/fine_pixel_asec)

target_coarse_res_asec = 6.5*60
coarse_bin_factor = int(np.round(target_coarse_res_asec/fine_pixel_asec))

# create analysis_root if doesn't exist
if not os.path.isdir(analysis_root):

	printbig('making directory'+analysis_root)
	os.mkdir(analysis_root)

# unpack .fits.fz to .fits if not already done
if not os.path.isfile(analysis_root + label + '_original.fits'):

	printbig('unpacking .fits into .fits.fz')
	print(raw_frames_root+label+'.fits.fz ==> ' + analysis_root+label+'_original.fits')
	call(['funpack','-E','1','-O',analysis_root+label+'_original.fits',raw_frames_root+label+'.fits.fz'])

# run swarp to grid the fine res stereographic frame to fine res orthographic
if not os.path.isfile(analysis_root+label+'.fits') or overwrite_swarp:

	printbig("swarping to %1.2f\" res ortho centered at (RA,Dec)=(%1.2f,%1.2f)"%(fine_pixel_asec,ra_cent_deg,dec_cent_deg))
	call(['swarp',analysis_root+label+'_original.fits','-c','grid_to_ortho.swarp', \
             '-IMAGEOUT_NAME',analysis_root+label+'.fits', \
             '-WEIGHTOUT_NAME','/dev/null', \
             '-RESAMPLING_TYPE','NEAREST', '-PIXEL_SCALE',str(fine_pixel_asec), \
             '-IMAGE_SIZE',str(n_fine)+','+str(n_fine),
             '-CENTER',str(ra_cent_deg)+','+str(dec_cent_deg)])

# generate source lists for masking
if not (os.path.isfile(analysis_root+label+'_artifacts_sources.reg') and \
           os.path.isfile(analysis_root+label+'_artifacts_sources.csv') and \
           os.path.isfile(analysis_root+label+'_artifacts.csv')) or overwrite_mask_source_lists:

	printbig('generating source lists for masking')
	call(['./catalog_to_mask_lists.py',analysis_root+label+'.fits',\
              raw_frames_root+'dph/'+label+'.dph',analysis_root,label])

# show masked regions on fits image in ds9 for debugging
if show_fits_image_with_masked_regions:
	call(['ds9',analysis_root+label+'.fits','-scale','zscale','-xpa','local','-regions',analysis_root+label+'_artifacts_sources.reg'])

# make fine res fits image of masks for debugging
# and with just bright pixel artifacts masked
if not (os.path.isfile(analysis_root+label+'_artifacts_sources.fits') and \
          os.path.isfile(analysis_root+label+'_artifacts.fits')) or overwrite_mask_fits_images:

	printbig('making fits images of masks')
	call(['./maskcsv2fits.py',analysis_root+label+'.fits', \
              analysis_root+label+'_artifacts_sources.csv', \
              analysis_root+label+'_artifacts_sources.fits'])	

        call(['./maskcsv2fits.py',analysis_root+label+'.fits', \
              analysis_root+label+'_artifacts.csv', \
              analysis_root+label+'_artifacts.fits'])

# bin to 3' resolution, excluding masked regions
if not (os.path.isfile(analysis_root+label+'_mask_artifacts_coarse.fits') and \
        os.path.isfile(analysis_root+label+'_mask_artifacts_sources_coarse.fits')) or override_coarse_binning:

	printbig('coarse binning by factor of '+str(coarse_bin_factor)+' to '+str(fine_pixel_asec*coarse_bin_factor/60.)+'\' res')
	call(['./coarse_bin_fits_image_with_mask.py',analysis_root+label+'.fits', \
             analysis_root+label+'_artifacts_sources.fits', \
             analysis_root+label+'_mask_artifacts_sources_coarse.fits', \
             analysis_root+label+'_mask_artifacts_sources_coarse_weights.fits',str(coarse_bin_factor)])

        call(['./coarse_bin_fits_image_with_mask.py',analysis_root+label+'.fits', \
             analysis_root+label+'_artifacts.fits', \
             analysis_root+label+'_mask_artifacts_coarse.fits', \
             analysis_root+label+'_mask_artifacts_coarse_weights.fits',str(coarse_bin_factor)])
	

