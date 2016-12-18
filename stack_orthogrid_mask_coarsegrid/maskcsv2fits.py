#!/Users/abrahamn/anaconda2/bin/python

import numpy as np
import commands
import sys
from astropy.io import fits
from astropy import wcs

# def mask_circle_in_image(img,xcent,ycent,r):
# 	n=img.shape[0]
# 	for y in range(-int(r),int(r)+1): # x is the x distance from circle center
# 		img[ycent+y,int(xcent-np.sqrt(r**2-y**2)):int(xcent+np.sqrt(r**2-y**2))] = 0

def mask_circle_in_image(img,xcent,ycent,r):
	n=img.shape[0]
	if not ((0<=xcent<n) and (0<=ycent<n)): return
	for y in range(-int(r),int(r)+1): # x is the x distance from circle center
		if not (0 <= ycent+y < n): continue
		#print(n)
		#print(int(xcent-np.sqrt(r**2-y**2)),int(xcent+np.sqrt(r**2-y**2)))
		#rint(max(0,int(xcent-np.sqrt(r**2-y**2))),min(n,int(xcent+np.sqrt(r**2-y**2))),)
		img[ycent+y,max(0,int(xcent-np.sqrt(r**2-y**2))):min(n,int(xcent+np.sqrt(r**2-y**2)))] = 0

if len(sys.argv) != 4:
	print('Usage: ./maskcsv2fits.py [original fits image] [CSV region file] [mask output fits file]')
	sys.exit(0)

original_img_path = sys.argv[1]
csvpath = sys.argv[2]
mask_fits_path = sys.argv[3]

regions = np.genfromtxt(csvpath,delimiter=',')
num_regions = len(regions)

commands.getoutput('cp '+original_img_path+' '+mask_fits_path)
hdulist = fits.open(mask_fits_path,mode='update')
h = hdulist[0].header
asec_per_pixel = np.abs(h['CD1_1'])*3600

w = wcs.WCS(h)
img = hdulist[0].data
img[:,:] = 1

# convert RA,Dec to pixel coordinates
px,py = w.wcs_world2pix(regions[:,0],regions[:,1],1)

# draw a filled circle of zeros on the mask image for each source
for i in range(num_regions):
	if i % 2500 == 0: print(1.*i/num_regions)
	mask_circle_in_image(img,px[i],py[i],regions[i,2]/asec_per_pixel)

hdulist.flush()
hdulist.close()
print('wrote '+mask_fits_path)


