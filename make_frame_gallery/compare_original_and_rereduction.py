import commands
import pyfits
import matplotlib.pyplot as plt

def drawbox(ra,dec,fov,col,lw):
	if ra > 180: ra -= 360
	plt.plot([ra-fov/2,ra-fov/2,ra+fov/2,ra+fov/2,ra-fov/2],[dec-fov/2,dec+fov/2,dec+fov/2,dec-fov/2,dec-fov/2],col+'-',linewidth=lw)

def generate_frame_context_image(rootdir,ra,dec,label):
	plt.figure(figsize=(3,3))
	drawbox(ra,dec,5.,'r',lw=2)
	drawbox(0.,-27.,10.,'k',lw=3)

	plt.xlim([-12.5,12.5])
	plt.ylim([-42.5,-17.5])

	plt.axes().set_aspect('equal', 'datalim')
	plt.savefig(rootdir+label+'_survey_context.png',bbox_inches='tight')

gal1 = '/Volumes/abraham/xcor_data/analysis/ATLAS_mwa57694/gallery/'
gal2 = '/Volumes/abraham/xcor_data/analysis/ATLAS_mwa57694_rereduction/gallery/'

imgnames = [l.split('/')[-1] for l in commands.getoutput('ls '+gal1+'*.jpg').split()]

f = open(gal2+'index.html','w')
f.write('<html><body><table style="text-align:center">\n')
f.write('<tr><td> </td><td><h2>ATLAS_mwa57694</h2></td><td><h2>ATLAS_mwa57694_rereduction</h2></td></tr>\n')

for imgname in imgnames:

	label = imgname.split('.')[0]
	# print(label)
	# hdulist = pyfits.open('/Volumes/abraham/xcor_data/ATLAS_mwa57694/' + label+'.fits')
	# hd = hdulist[0].header
	# ra = hd['RA-MNT']
	# dec = hd['DEC-MNT']

	#generate_frame_context_image(gal2,ra,dec,label)

	f.write('<tr>\n')
	f.write('<td>'+label+'<br/><img src="'+label+'_survey_context.png"/>\n</td>')
	f.write('<td><img src="ATLAS_mwa57694/'+imgname+'"/></td>\n')
	f.write('<td><img src="ATLAS_mwa57694_rereduction/'+imgname+'"/></td>\n')
	f.write('</tr>\n')


f.write('</table></body></html>')
f.close()


