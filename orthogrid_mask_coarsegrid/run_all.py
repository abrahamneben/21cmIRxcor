#!/Users/abrahamn/anaconda2/bin/python

from subprocess import call

labels = ['02a57639o0352I','02a57639o0356I','02a57639o0342I','02a57639o0346I']

for label in labels:
	call(['./run.py',label])

