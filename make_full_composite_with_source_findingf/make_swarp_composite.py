#!/Users/abrahamn/anaconda2/bin/python

import commands
from subprocess import call


# find all labels
raw_frames_root = '/Volumes/abraham/xcor_data/ATLAS_mwa57694_rereduction/'
labels = [l.split('/')[-1].split('.')[0] for l in commands.getoutput('ls '+raw_frames_root+'*.fits.fz').split()]

# funpack everything
# for l in labels:
# 	print(raw_frames_root+l+'.fits.fz ==> ' + raw_frames_root+l+'.fits')
# 	call(['funpack','-E','1','-O',raw_frames_root+l+'.fits',raw_frames_root+l+'.fits.fz'])

# run swarp
analysis_root = '/Volumes/abraham/xcor_data/analysis/ATLAS_mwa57694_rereduction/full_composite/'
analysis_name = 'full_composite'

ra_cent_deg = 0.
dec_cent_deg = -27.
fine_pixel_asec = 1.86
fov_deg = 15
n_fine = int(fov_deg*3600./fine_pixel_asec)

call(['swarp']+[raw_frames_root+l+'.fits' for l in labels]+['-c','full_composite.swarp', \
'-IMAGEOUT_NAME',analysis_root+analysis_name+'.fits', \
'-WEIGHTOUT_NAME',analysis_root+analysis_name+'_weight.fits', \
'-RESAMPLING_TYPE','NEAREST', '-PIXEL_SCALE',str(fine_pixel_asec), \
'-IMAGE_SIZE',str(n_fine)+','+str(n_fine),
'-CENTER',str(ra_cent_deg)+','+str(dec_cent_deg)])

# run sextractor
# call(['sex',analysis_root+analysis_name+'.fits','-c','full_composite.sex'])


