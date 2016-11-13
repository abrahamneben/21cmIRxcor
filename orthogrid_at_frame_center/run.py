#!/Users/abrahamn/anaconda2/bin/python

import os
from subprocess import call
import numpy as np
import sys
from astropy.io import fits

def printbig(s):
	n = len(s)
	print('')
	print('*'*(n+4))
	print('* '+s+' *')
	print('*'*(n+4))
	print('')


if len(sys.argv) != 3:
	print('Usage: ./run.py [ATLAS_mwa57694_rereduction] [02a57639o0338I]')
	sys.exit(0)

run = sys.argv[1]
label = sys.argv[2]
raw_frames_root = '/Volumes/abraham/xcor_data/'+run+'/'
analysis_root = '/Volumes/abraham/xcor_data/analysis/' + run + '/' + label + '/'
printbig('processing '+label)

hdulist = fits.open(raw_frames_root+label+'.fits')
h = hdulist[0].header
ra_cent_deg = h['RA']
dec_cent_deg = h['DEC']

fine_pixel_asec = 1.86 # same as original stereographic projection res
fov_deg = 6.
n_fine = int(fov_deg*3600/fine_pixel_asec)

if not os.path.exists(analysis_root):
  os.makedirs(analysis_root)

# run swarp to grid the fine res stereographic frame to fine res orthographic
printbig("swarping to %1.2f\" res ortho centered at (RA,Dec)=(%1.2f,%1.2f)"%(fine_pixel_asec,ra_cent_deg,dec_cent_deg))
call(['swarp',raw_frames_root+label+'.fits','-c','grid_to_ortho.swarp', \
             '-IMAGEOUT_NAME',analysis_root+label+'_5degframecentered.fits', \
             '-WEIGHTOUT_NAME','/dev/null', \
             '-RESAMPLING_TYPE','NEAREST', '-PIXEL_SCALE',str(fine_pixel_asec), \
             '-IMAGE_SIZE',str(n_fine)+','+str(n_fine),
             '-CENTER',str(ra_cent_deg)+','+str(dec_cent_deg)])

