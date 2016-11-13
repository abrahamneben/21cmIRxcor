#!/Users/abrahamn/anaconda2/bin/python

from subprocess import call

runs = ['ATLAS_mwa57694','ATLAS_mwa57694_rereduction']
labels = ['02a57694o0313I','02a57694o0313I']

for i in range(len(labels)):
	call(['./run.py',runs[i],labels[i]])

