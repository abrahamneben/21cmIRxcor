
import numpy as np
import commands

def obsids2txt(obsids,fname):
	f=open(fname,'w')
	f.write('\n'.join(map(str,list(obsids)) ))
	f.close()

def get_jd_from_obsid(o):
	return commands.getoutput('get_observation_info.py -g '+str(1069761200)).split()[0].split('_')[-1]

wedgecut_and_rescut_obsids_path = '/nfs/eor-00/h1/beards/obs_lists/long_runs/wedge_cut_plus_res_cut.txt'
obsids = np.genfromtxt(wedgecut_and_rescut_obsids_path,dtype=int)
num_obsids = len(obsids)

jds = [get_jd_from_obsid(o) for o in obsids]

#obsids2txt(obsids_half_1,'obsids_wedgecutrescut_half_1.txt')


