#!/Users/abrahamn/anaconda2/bin/python

from subprocess import call

analysis_root = '/Volumes/abraham/xcor_data/analysis/ATLAS_mwa57639/whole_field/'
analysis_name = 'whole_field_bgsub'

# run sextractor
call(['sex',analysis_root+analysis_name+'.fits','-c','source_finding.sex',\
		  '-CATALOG_NAME',analysis_root+analysis_name+'.sex',])
