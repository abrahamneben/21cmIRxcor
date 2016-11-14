#!/Users/abrahamn/anaconda2/bin/python

from subprocess import call

for fieldi in range(4):

	# run swarp
	raw_frames_root = '/Volumes/abraham/xcor_data/ATLAS_mwa57694_rereduction/'
	analysis_root = '/Volumes/abraham/xcor_data/analysis/ATLAS_mwa57694_rereduction/field'+str(fieldi)+'/'
	analysis_name = 'field'+str(fieldi)

	# run sextractor
	call(['sex',analysis_root+analysis_name+'.fits','-c','source_finding.sex',\
		  '-CATALOG_NAME',analysis_root+analysis_name+'.sex',])

	# run sextractor
# call(['sex','/Volumes/abraham/xcor_data/ATLAS_mwa57694_rereduction/02a57694o0275I.fits',\
# 	'-c','source_finding.sex','-CATALOG_NAME','test.sextractor'])

