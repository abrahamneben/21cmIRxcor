#unpack .fits.fz to .fits if not already done

import commands
import sys
from subprocess import call

raw_frames_root = '/volumes/abraham/xcor_data/ATLAS_mwa57639/'
labels = [l.split('/')[-1].split('.')[0] for l in commands.getoutput('ls '+raw_frames_root+'*.fits.fz').split()]

print(labels)

for l in labels:
	print(raw_frames_root+l+'.fits.fz ==> ' + raw_frames_root+l+'.fits')
	call(['funpack','-E','1','-O',raw_frames_root+l+'.fits',raw_frames_root+l+'.fits.fz'])
