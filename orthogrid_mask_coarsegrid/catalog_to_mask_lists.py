#!/Users/abrahamn/anaconda2/bin/python

# given a dophot catalog and fits image, generate .reg and .csv
# source lists for masking. The lists contain sources from the dophot
# catalog, as well as the bad pixel artifacts (ie, the pixels that are way too bright)
# found in the fits image in pixel coordinates, then transformed using astroy to RA,Dec

import sys
import numpy as np
from astropy.io import fits
from astropy import wcs

if len(sys.argv) != 5:
	print('Usage: ./catalog_to_mask_lists.py [fits image] [dophot catalog] [output dir] [label]')
	sys.exit(0)

ortho_img_fits_path = sys.argv[1] # doesn't seem to work with the original stereographic image,
				  # only with the ortho image made by swarp
dophot_cat_path = sys.argv[2]
output_dir = sys.argv[3]
label = sys.argv[4]

source_mask_radius_asec = 30.
bad_pixel_mask_radius_asec = 150.

bad_pixel_threshold_adu = 17.e3

# add a slash to the end of output_dir if n eeded
if output_dir[-1] != '/': output_dir += '/'

# open files for writing
artifacts_sources_reg_file = open(output_dir + label + '_artifacts_sources.reg','w')
artifacts_sources_csv_file = open(output_dir + label + '_artifacts_sources.csv','w')
artifacts_csv_file = open(output_dir + label + '_artifacts.csv','w')

artifacts_sources_reg_file.write('# Region file format: DS9 version 4.1\n')
artifacts_sources_reg_file.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
artifacts_sources_reg_file.write('icrs\n')

# mask a fixed radius arounda all dophot sources
srcdat = np.genfromtxt(dophot_cat_path)
for i in range(len(srcdat)):
	artifacts_sources_reg_file.write("circle(%f,%f,%f\")\n"%(srcdat[i,0],srcdat[i,1],source_mask_radius_asec))
	artifacts_sources_csv_file.write("%f,%f,%f\n"%(srcdat[i,0],srcdat[i,1],source_mask_radius_asec))


# now look for bad (ie, nearly saturated) pixels and mask big circles around them to remove artifacts
print('opening '+ortho_img_fits_path+'...')
hdulist = fits.open(ortho_img_fits_path)
img = hdulist[0].data
w = wcs.WCS(hdulist[0].header)
hdulist.close()

xinds,yinds = np.where(img>bad_pixel_threshold_adu)
radec = w.wcs_pix2world(np.array([yinds,xinds]).T, 0)

for i in range(len(xinds)):
	artifacts_sources_reg_file.write("circle(%f,%f,%f\") # color=blue\n"%(radec[i,0],radec[i,1],bad_pixel_mask_radius_asec))
	artifacts_sources_csv_file.write("%f,%f,%f\n"%(radec[i,0],radec[i,1],bad_pixel_mask_radius_asec))
	artifacts_csv_file.write("%f,%f,%f\n"%(radec[i,0],radec[i,1],bad_pixel_mask_radius_asec))

# close and save ds9 .reg file
artifacts_sources_reg_file.close()
artifacts_sources_csv_file.close()
artifacts_csv_file.close()

print('wrote '+output_dir+label+'_artifacts_sources.reg')
print('wrote '+output_dir+label+'_artifacts_sources.csv')
print('wrote '+output_dir+label+'_artifacts.csv')

