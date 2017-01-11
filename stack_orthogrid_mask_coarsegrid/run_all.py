#!/Users/abrahamn/anaconda2/bin/python

from subprocess import call

#runs = ['ATLAS_mwa57694','ATLAS_mwa57694_rereduction']
runs = ['ATLAS_mwa57694_rereduction']

label_groups = [ \
'02a57694o0299I 02a57694o0303I 02a57694o0307I 02a57694o0319I 02a57694o0323I 02a57694o0327I 02a57694o0339I 02a57694o0343I 02a57694o0347I', \
'02a57694o0298I 02a57694o0302I 02a57694o0306I 02a57694o0318I 02a57694o0322I 02a57694o0326I 02a57694o0338I 02a57694o0342I 02a57694o0346I',\
'02a57694o0301I 02a57694o0305I 02a57694o0309I 02a57694o0321I 02a57694o0325I 02a57694o0329I 02a57694o0341I 02a57694o0345I 02a57694o0349I',\
'02a57694o0300I 02a57694o0304I 02a57694o0308I 02a57694o0320I 02a57694o0324I 02a57694o0328I 02a57694o0340I 02a57694o0344I 02a57694o0348I',\
# '02a57639o0347I 02a57639o0348I 02a57639o0349I 02a57639o0350I 02a57639o0351I', \
# '02a57639o0352I 02a57639o0353I 02a57639o0354I 02a57639o0355I 02a57639o0356I', \
# '02a57639o0337I 02a57639o0338I 02a57639o0339I 02a57639o0340I 02a57639o0341I' \
]

for run in runs:
	for label_group in label_groups:
		call(['caffeinate','-i','./run.py',run]+label_group.split(' '))

