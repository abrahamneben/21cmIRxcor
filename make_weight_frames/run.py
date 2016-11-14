#!/Users/abrahamn/anaconda2/bin/python

import commands
from astropy.io import fits
import numpy as np

raw_frames_root = '/Volumes/abraham/xcor_data/ATLAS_mwa57694_rereduction/'
labels = [l.split('/')[-1].split('.')[0] for l in commands.getoutput('ls '+raw_frames_root+'*.fits.fz').split()]

for l in labels:
	print(l)
	commands.getoutput('cp '+raw_frames_root+l+'.fits '+raw_frames_root+l+'_weight.fits')

	hdulist = fits.open(raw_frames_root+l+'_weight.fits',mode='update')
	hdulist[0].data = np.uint16(hdulist[0].data > 200)

	hdulist.flush()
	hdulist.close()

