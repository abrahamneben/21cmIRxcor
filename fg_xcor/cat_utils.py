import numpy as np
from astropy.io import fits
from astropy import wcs
from scipy.io import readsav

class Catalog:
    def __init__(self,ra,dec,jy):
        self.ra = ra
        self.dec = dec
        self.jy = jy
         
    def calc_min_max_ra_dec(self):
        self.min_ra, self.max_ra, self.min_dec, self.max_dec = np.min(self.ra),np.max(self.ra),np.min(self.dec),np.max(self.dec)
        
    def limit_to_ra_dec_min_max_of_other_cat(self,cat):
        g1 = (self.ra > cat.min_ra)&(self.ra < cat.max_ra)&(self.dec > cat.min_dec)&(self.dec < cat.max_dec)
        return Catalog(self.ra[g1],self.dec[g1],self.jy[g1])

class MWACatalog(Catalog):
    def __init__(self,fhdsavpath):
        self.jy2pointsize = 10.
        self.col = 'r'
        
        mwadat = readsav(fhdsavpath)
        # mwadat['catalog'].dtype.names

        self.ra, self.dec = mwadat['catalog'].ra,mwadat['catalog'].dec
        self.ra[self.ra>180] -= 360
        self.calc_min_max_ra_dec()
        
        self.jy = np.array([f.i[0] for f in mwadat['catalog'].flux])
    
class IRCatalog(Catalog):
    def __init__(self,raw_frames_path,label,orthofitsimagepath):
        print('loading '+label)
        self.jy2pointsize = 500.
        self.col = 'b'
        
        irdat = np.genfromtxt(raw_frames_path+label+'.dph',usecols=(0,1,2))

        self.ra_all = irdat[:,0]
        self.ra_all[self.ra_all>180] -= 360
        self.dec_all = irdat[:,1]
        
        self.m_all = irdat[:,2]
        self.jy_all = 3631.*10.**(-self.m_all/2.5)
        
        self.identify_catalog_artifacts(orthofitsimagepath)
        self.calc_min_max_ra_dec()
        
    def identify_catalog_artifacts(self,orthofitsimagepath):
        print('identifying and excluding artifacts (ie, saturated pixels)')
        bad_pixel_threshold_adu = 17.e3
        artifact_masking_radius_asec = 80.
        hdulist = fits.open(orthofitsimagepath)

        img = hdulist[0].data
        w = wcs.WCS(hdulist[0].header)
        hdulist.close()

        xinds,yinds = np.where(img>bad_pixel_threshold_adu)
        radec = w.wcs_pix2world(np.array([yinds,xinds]).T, 0)
        art_ra = radec[:,0]
        art_ra[art_ra>180] -= 360
        art_dec = radec[:,1]

        num_artifacts_near_every_source = np.zeros(len(self.ra_all),dtype=bool) # equals true if the source is not an artifact
        for i in range(len(art_ra)):
            num_artifacts_near_every_source += np.sqrt((self.dec_all-art_dec[i])**2+np.cos(self.dec_all)**2*(self.ra_all-art_ra[i])**2) < artifact_masking_radius_asec/3600
        g = num_artifacts_near_every_source == 0
        
        self.ra = self.ra_all[g]
        self.dec = self.dec_all[g]
        self.jy = self.jy_all[g]
            
def plot_cat_list(plt,cats,jymin=.2,alpha=.2,exclude_artifacts=True):
    for cat in cats:
        ra,dec,jy = (cat.ra,cat.dec,cat.jy) if exclude_artifacts else (cat.ra_all,cat.dec_all,cat.jy_all)
        plt.scatter(ra[jy>jymin],\
                    dec[jy>jymin],\
                    cat.jy2pointsize*jy[jy>jymin],cat.col,alpha=alpha)
        plot_catalog_boundary(plt,cat)
        
def plot_catalog_boundary(plt,cat):
    ra0,ra1,dec0,dec1 = cat.min_ra,cat.max_ra,cat.min_dec,cat.max_dec
    plt.plot([ra0,ra1,ra1,ra0,ra0],[dec0,dec0,dec1,dec1,dec0],'-')
    
    
    
    