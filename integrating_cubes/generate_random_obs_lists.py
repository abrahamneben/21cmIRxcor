
import numpy as np
from numpy.random import choice

def obsids2txt(obsids,fname):
	f=open(fname,'w')
	f.write('\n'.join(map(str,list(obsids)) ))
	f.close()

wedgecut_and_rescut_obsids_path = '/nfs/eor-00/h1/beards/obs_lists/long_runs/wedge_cut_plus_res_cut.txt'
obsids = np.genfromtxt(wedgecut_and_rescut_obsids_path,dtype=int)
num_obsids = len(obsids)

obsids_half_1 = choice(obsids,size=num_obsids/2,replace=False)
obsids_half_2 = np.array(list(set(obsids)-set(obsids_half_1)))

obsids_quarter_1 = choice(obsids,size=num_obsids/4,replace=False)
obsids_quarter_2 = choice(np.array(list(set(obsids)-set(obsids_quarter_1))),size=num_obsids/4,replace=False)
obsids_quarter_3 = choice(np.array(list(set(obsids)-set(obsids_quarter_1)-set(obsids_quarter_2))),size=num_obsids/4,replace=False)
obsids_quarter_4 = np.array(list(set(obsids)-set(obsids_quarter_1)-set(obsids_quarter_2)-set(obsids_quarter_3)))

obsids2txt(obsids_half_1,'obsids_wedgecutrescut_half_1.txt')
obsids2txt(obsids_half_2,'obsids_wedgecutrescut_half_2.txt')

obsids2txt(obsids_quarter_1,'obsids_wedgecutrescut_quarter_1.txt')
obsids2txt(obsids_quarter_2,'obsids_wedgecutrescut_quarter_2.txt')
obsids2txt(obsids_quarter_3,'obsids_wedgecutrescut_quarter_3.txt')
obsids2txt(obsids_quarter_4,'obsids_wedgecutrescut_quarter_4.txt')


