#!/Users/abrahamn/anaconda2/bin/python

import sys
from astropy.io import fits
import os

# from http://scipy.github.io/old-wiki/pages/Cookbook/Rebinning, Example 2
def bin_by_factor(img,n_coarse,factor):
	return img.reshape(n_coarse,factor,n_coarse,factor,).sum(1).sum(2)

if len(sys.argv) != 6:
	print('Usage: ./coarse_bin_fits_image_with_mask.py [fine res fits image] [fine res fits mask image] [output fits image] [output fits weights] [coarse bin factor]')
	sys.exit(0)

img_fits_path = sys.argv[1]
mask_fits_path = sys.argv[2]
output_fits_path = sys.argv[3]
output_fits_weights_path = sys.argv[4]
coarse_bin_factor = int(sys.argv[5])

# read in img
hdulist_img = fits.open(img_fits_path)
img = hdulist_img[0].data
h = hdulist_img[0].header
n = img.shape[0]
hdulist_img.close()

# read in mask
hdulist_mask = fits.open(mask_fits_path)
mask = hdulist_mask[0].data
hdulist_mask.close()

n_coarse = int(n/coarse_bin_factor)
n_fine_max = n_coarse*coarse_bin_factor

print('binning image')
masked_img_binned = bin_by_factor((img*mask)[:n_fine_max,:n_fine_max],n_coarse,coarse_bin_factor)
print('binning mask')
mask_binned = bin_by_factor(mask[:n_fine_max,:n_fine_max],n_coarse,coarse_bin_factor)
print('normalizing image')
masked_img_binned_normed = masked_img_binned / mask_binned * coarse_bin_factor**2
mask_binned[masked_img_binned == 0] = 0

# make a new fits header for the coarse image
h2 = h
h2['CD1_1'] *= coarse_bin_factor
h2['CD2_2'] *= coarse_bin_factor
h2['NAXIS1'] = n_coarse
h2['NAXIS2'] = n_coarse
h2['CRPIX1'] = n_coarse/2
h2['CRPIX2'] = n_coarse/2

# remove old files is they exist
if os.path.isfile(output_fits_path): os.remove(output_fits_path)
if os.path.isfile(output_fits_weights_path): os.remove(output_fits_weights_path)

# write the coarse gridded image to new fits file
print('writing '+output_fits_path)
hdu_coarse_img = fits.PrimaryHDU(masked_img_binned_normed)
hdu_coarse_img.header = h2
hdu_coarse_img.writeto(output_fits_path,output_verify='warn')

# write coarse gridded weights to new fits file
print('writing '+output_fits_weights_path)
hdu_coarse_weights = fits.PrimaryHDU(mask_binned)
hdu_coarse_weights.header = h2
hdu_coarse_weights.writeto(output_fits_weights_path, output_verify='warn')

