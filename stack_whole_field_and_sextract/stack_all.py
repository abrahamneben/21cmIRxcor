from subprocess import call

raw_frames_root = '/volumes/abraham/xcor_data/ATLAS_mwa57639/'
analysis_root = '/volumes/abraham/xcor_data/analysis/ATLAS_mwa57639/'
analysis_name = 'whole_field'

labels = open('ATLAS_mwa57639_labels_good.txt').read().split('\n')

print(labels)

fine_pixel_asec = 1.86
fov_deg = 5.
n_fine = int(fov_deg*3600./fine_pixel_asec)

ra_cent_deg = 0.
dec_cent_deg = -30.

fine_pixel_asec = 1.86
fov_deg = 20.
n_fine = int(fov_deg*3600./fine_pixel_asec)

print(n_fine)

call(['swarp']+[raw_frames_root+l+'.fits' for l in labels]+\
['-WEIGHT_TYPE','MAP_WEIGHT']+\
['-WEIGHT_IMAGE']+[' '.join([raw_frames_root+l+'_weight.fits' for l in labels])]+\
['-c','stack_whole_field.swarp', \
'-IMAGEOUT_NAME',analysis_root+analysis_name+'_bgsub.fits', \
'-WEIGHTOUT_NAME',analysis_root+analysis_name+'_bgsub_weight.fits', \
'-RESAMPLING_TYPE','NEAREST', '-PIXEL_SCALE',str(fine_pixel_asec), \
'-IMAGE_SIZE',str(n_fine)+','+str(n_fine),
'-CENTER',str(ra_cent_deg)+','+str(dec_cent_deg)])

