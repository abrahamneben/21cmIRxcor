
import numpy as np
from astropy.time import Time
import commands
import pyfits

import matplotlib.pyplot as plt

col1 = '#fff'
col2 = '#ddd'
def flipcolor(currcolor):
	if currcolor == col1: return col2
	return col1

def drawbox(ra,dec,fov,col,lw):
	if ra > 180: ra -= 360
	plt.plot([ra-fov/2,ra-fov/2,ra+fov/2,ra+fov/2,ra-fov/2],[dec-fov/2,dec+fov/2,dec+fov/2,dec-fov/2,dec-fov/2],col+'-',linewidth=lw)

def generate_frame_context_image(ra,dec,label):
	plt.figure(figsize=(3,3))
	drawbox(ra,dec,5.,'r',lw=2)
	drawbox(0.,-27.,15.,'k',lw=3)

	plt.xlim([-12.5,12.5])
	plt.ylim([-42.5,-17.5])

	plt.axes().set_aspect('equal', 'datalim')
	plt.savefig('images/'+label+'_survey_context.png',bbox_inches='tight')

def plot_all_frames_context_image(ra_centers,dec_centers):
	ra_centers[ra_centers>180] -= 360

	plt.figure(figsize=(3,3))

	for i in range(len(ra_centers)):
		drawbox(ra_centers[i],dec_centers[i],5.,'r',lw=2)

	drawbox(0.,-27,15.,'k',lw=3)

	plt.savefig('images/all_frames_context.png',bbox_inches='tight')

raw_data_root = '/Volumes/abraham/xcor_data/ATLAS_mwa57639/'
atlas_fnames = commands.getoutput('ls images/*.jpg').split('\n')

f = open('index.html','w')
f.write('<html><head>')
f.write('<style>p{margin-bottom:0;margin-top:5px}</style>')
f.write('<title>ATLAS 57639, MWA EOR0 field</title></head><body>')

f.write('<p>all images are shown on a scale of 112*skyadu to 118*skyadu</p>')
f.write('<table width="100%">')

currcolor = col1
ra_old = 0
dec_old = 0

ra_centers = []
dec_centers = []

f.write('<tr>\n')
for i in range(len(atlas_fnames)):
	label = atlas_fnames[i].split('/')[1].split('_')[0]
	print(label)

	hdulist = pyfits.open(raw_data_root + label+'.fits.fz')
	hd = hdulist[1].header
	
	ra = hd['RA-MNT']
	dec = hd['DEC-MNT']
	az = str(hd['AZ'])
	alt = str(hd['ALT'])
	dt = str(Time(hd['MJD-OBS'], format='mjd').iso)
	nstar = str(hd['NSTAR'])
	sky_frac_rms = str(1.*hd['SKYRMS']/hd['SKYADU'])
	magzpt = str(hd['MAGZPT'])
	hdulist.close()		

	dtheta = np.sqrt((ra-ra_old)**2+(dec-dec_old)**2)
	
	if dtheta>.5:
		currcolor = flipcolor(currcolor)
		f.write('</tr>\n<tr style="background-color:'+currcolor+'">\n')

		ra_centers.append(ra)
		dec_centers.append(dec)

		#generate_frame_context_image(ra,dec,label)
		f.write('<td><img src="images/'+label+'_survey_context.png"/>\n</td>')

	f.write('<td>\n\t<p><b>'+label+'</b></p>\n\t<p>(RA,Dec)=('+str(ra)+','+str(dec)+')</p>\n\t<p>(Alt,Az)=('+alt+','+az+')</p>\n\t<p>'+dt+'</p>')
	f.write('<p>nstar in tphot: '+nstar+'</p><p>skyrms/skyadu = '+sky_frac_rms+'</p><p>magzpt = '+magzpt+'</p>')
	f.write('<p><img width="100%" float="left" src="images/'+label+'_20asec.jpg"/></p>\n</td>\n')

	ra_old = ra
	dec_old = dec

f.write('</tr>\n</table>')

plot_all_frames_context_image(ra_centers,dec_centers)
f.write('<p><img src="images/all_frames_context.png"/></p>')

f.write('<p><img src="images/0372-0376.gif"/></p>')

f.write('</body></html>')
f.close()

