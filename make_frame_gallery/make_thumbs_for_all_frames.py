
import os
import commands
from astropy.io import fits
from subprocess import call

def printbig(s):
        n = len(s)
        print('')
        print('*'*(n+4))
        print('* '+s+' *')
        print('*'*(n+4))
        print('')

raw_frames_root = '/Volumes/abraham/xcor_data/ATLAS_mwa57694_rereduction/'
analysis_root = '/Volumes/abraham/xcor_data/analysis/ATLAS_mwa57694_rereduction/gallery/'


# create analysis_root if doesn't exist
if not os.path.isdir(analysis_root):
        printbig('making directory'+analysis_root)
        os.mkdir(analysis_root)

# find all labels
labels = [l.split('/')[-1].split('.')[0] for l in commands.getoutput('ls '+raw_frames_root+'*.fits.fz').split()]
#print(labels)

# unpack 
for l in labels:
	print('unpacking '+l)
	call(['funpack','-E','1','-O',raw_frames_root+l+'.fits',raw_frames_root+l+'.fits.fz'])


# # make thumbnails
for l in labels:
	print(l)

	# find sky background level
	hdulist = fits.open(raw_frames_root+l+'.fits')
	h = hdulist[0].header
	skyadu = h['SKYADU']

	call(['ds9',raw_frames_root+l+'.fits','-scale','limits',\
		str(.9*skyadu),str(1.1*skyadu),'-width','600','-xpa','local','-zoom','to','fit',\
		'-saveimage',analysis_root+l+'.jpg','-exit'])
