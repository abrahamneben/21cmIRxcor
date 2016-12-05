import numpy as np
from astropy.io import fits
from astropy import wcs
from scipy.io import readsav
import scipy.signal

class Catalog:
    def __init__(self,ra,dec,jy):
        self.ra = ra
        self.dec = dec
        self.jy = jy
        self.calc_min_max_ra_dec()
         
    def calc_min_max_ra_dec(self):
        self.min_ra, self.max_ra, self.min_dec, self.max_dec = np.min(self.ra),np.max(self.ra),np.min(self.dec),np.max(self.dec)
        self.mean_dec = np.mean(self.dec)
        
    def limit_to_ra_dec_min_max_of_other_cat(self,cat):
        g1 = (self.ra > cat.min_ra)&(self.ra < cat.max_ra)&\
             (self.dec > cat.min_dec)&(self.dec < cat.max_dec)
        return Catalog(self.ra[g1],self.dec[g1],self.jy[g1])
    
    def limit_to_ra_dec_min_max(self,min_ra,max_ra,min_dec,max_dec):
        g1 = (self.ra > min_ra)&(self.ra < max_ra)&\
             (self.dec > min_dec)&(self.dec < max_dec)
        return Catalog(self.ra[g1],self.dec[g1],self.jy[g1])
    
    def join_with_cat(self,cat):
        return Catalog(np.append(self.ra,cat.ra), np.append(self.dec,cat.dec), np.append(self.jy,cat.jy))
    

class MWACatalog(Catalog):
    def __init__(self,fhdsavpath):
        mwadat = readsav(fhdsavpath)
        # mwadat['catalog'].dtype.names

        self.ra, self.dec = mwadat['catalog'].ra,mwadat['catalog'].dec
        self.ra[self.ra>180] -= 360
        self.calc_min_max_ra_dec()
        
        self.jy = np.array([f.i[0] for f in mwadat['catalog'].flux])
    
class IRCatalog(Catalog):
    # note that the given fits image must be in ortho projection, ie, with swarp
    def __init__(self,dph_path='',se_path='',fits_path='',se_magzpt=20.56):
        assert dph_path != '' or se_path != ''
        
        if dph_path != '':
            print('loading '+dph_path)
            irdat = np.genfromtxt(dph_path,usecols=(0,1,2))

            self.ra_all = irdat[:,0]
            self.ra_all[self.ra_all>180] -= 360
            self.dec_all = irdat[:,1]

            self.m_all = irdat[:,2]
            self.jy_all = 3631.*10.**(-self.m_all/2.5)
        elif se_path != '':
            print('loading'+se_path)
            se_dat = np.genfromtxt(se_path)
            se_counts = se_dat[:,1]
            self.ra_all = se_dat[:,3]
            self.ra_all[self.ra_all>180] -= 360
            self.dec_all = se_dat[:,4]
            self.jy_all = 3631*10.**(-se_magzpt/2.5)*se_counts/30
        
        if fits_path != '':
            self.identify_catalog_artifacts(fits_path)
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

        num_artifacts_near_every_source = np.zeros(len(self.ra_all),dtype=bool)
        for i in range(len(art_ra)):
            if i % 1000 == 0: print(1.*i/len(art_ra))
            num_artifacts_near_every_source += np.sqrt((self.dec_all-art_dec[i])**2+\
                                                       np.cos(self.dec_all)**2*(self.ra_all-art_ra[i])**2) < artifact_masking_radius_asec/3600
        g = num_artifacts_near_every_source == 0
        
        self.ra = self.ra_all[g]
        self.dec = self.dec_all[g]
        self.jy = self.jy_all[g]
            
def plot_cat_list(plt,cats,jy2pointsize_arr,cols,jymin_vals,alpha_vals,exclude_artifacts=True):
    for i in range(len(cats)):
        cat = cats[i]
        ra,dec,jy = (cat.ra,cat.dec,cat.jy) if exclude_artifacts else (cat.ra_all,cat.dec_all,cat.jy_all)
        plt.scatter(ra[jy>jymin_vals[i]],\
                    dec[jy>jymin_vals[i]],\
                    jy2pointsize_arr[i]*jy[jy>jymin_vals[i]],cols[i],alpha=alpha_vals[i])
        plot_catalog_boundary(plt,cat,cols[i])
        
def plot_catalog_boundary(plt,cat,col):
    ra0,ra1,dec0,dec1 = cat.min_ra,cat.max_ra,cat.min_dec,cat.max_dec
    plt.plot([ra0,ra1,ra1,ra0,ra0],[dec0,dec0,dec1,dec1,dec0],col+'-')
    
def cat2img(cat,bound_cat,dtheta,jymin=0,jymax=1.e9,verbose=False):
    fluxcut = (cat.jy < jymax)&(cat.jy > jymin)
    thetax_fluxcut = np.cos(cat.mean_dec*np.pi/180.)*cat.ra[fluxcut]
    thetay_fluxcut = cat.dec[fluxcut]
    jy_fluxcut = cat.jy[fluxcut]

    thetax0 = np.cos(bound_cat.mean_dec*np.pi/180.)*bound_cat.min_ra
    thetay0 = bound_cat.min_dec

    fov = np.min([bound_cat.max_dec-bound_cat.min_dec,np.cos(bound_cat.mean_dec*np.pi/180.)*(bound_cat.max_ra-bound_cat.min_ra)])
    n = int(fov/dtheta)

    img = np.zeros((n,n))
    for xi in range(n):
        if verbose and xi % 20 == 0: print(1.*xi/n)
        for yi in range(n):
            inpixel = (thetax_fluxcut > thetax0+xi*dtheta)&(thetax_fluxcut < thetax0+(xi+1)*dtheta)&\
                      (thetay_fluxcut > thetay0+yi*dtheta)&(thetay_fluxcut < thetay0+(yi+1)*dtheta)
            img[xi,yi] = np.sum(jy_fluxcut[inpixel])
    
    return img

def make_hann(n):
    w = scipy.signal.hann(n)
    wx,wy = np.meshgrid(w,w)
    w2 = wx*wy
    return w2, np.sqrt(np.mean(w2**2))

def calc_xspec(img1,img2,dtheta_deg,nbins,lmin,lmax,hann=True):
    assert img1.shape == img2.shape
    
    n,dang = img1.shape[0],dtheta_deg*np.pi/180.
    
    hann2D,hann2Drms = np.ones((n,n)),1
    if hann:
        hann2D,hann2Drms = make_hann(n)
    lvals = np.fft.fftfreq(n)*2*np.pi/dang
    lx,ly = np.meshgrid(lvals,lvals)
    lmag  = np.sqrt(lx**2+ly**2)

    img1_ft = np.fft.fft2(hann2D*(img1-img1.mean()))/hann2Drms
    img2_ft = np.fft.fft2(hann2D*(img2-img2.mean()))/hann2Drms
    
    lbinedges = np.linspace(lmin,lmax,nbins+1)
    lbincenters = .5*(lbinedges[0:nbins]+lbinedges[1:nbins+1])
    
    xspec_binned = np.zeros(nbins)
    pspec1_binned = np.zeros(nbins)
    pspec2_binned = np.zeros(nbins)
    bin_counts = np.zeros(nbins)
    for bini in range(nbins):
        inbin = (lmag>lbinedges[bini])&(lmag<lbinedges[bini+1])
        bin_counts[bini] = np.sum(inbin)

        xspec_binned[bini] = np.real(np.mean(img1_ft[inbin]*np.conj(img2_ft[inbin])))
        pspec1_binned[bini] = np.mean(np.abs(img1_ft[inbin])**2)
        pspec2_binned[bini] = np.mean(np.abs(img2_ft[inbin])**2)
    
    pspec_norm = (dang**2)/(n**2)
    return lbincenters,pspec1_binned*pspec_norm,pspec2_binned*pspec_norm,xspec_binned*pspec_norm,bin_counts

def plot_spectra(plt,lbins,pspec1,pspec2,xspec,bin_counts,ir_mwa_jymin_max,irlim=[1.e-11,1.e-8],mwalim=[1.e-12,1.e-6]):
    ir_jymin,ir_jymax,mwa_jymin,mwa_jymax = ir_mwa_jymin_max
    
    plt.figure(figsize=(20,4))
    
    plt.subplot(141)
    plt.loglog(lbins,pspec1)
    plt.xlim([np.min(lbins),np.max(lbins)])
    plt.ylim(irlim)
    plt.xlabel('\ell')
    plt.ylabel('C_\ell')
    plt.title("IR, (%1.0e < flux < %1.0e) jy"%(ir_jymin,ir_jymax))
    
    plt.subplot(142)
    plt.loglog(lbins,pspec2)
    plt.xlim([np.min(lbins),np.max(lbins)])
    plt.ylim(mwalim)
    plt.xlabel('\ell')
    plt.ylabel('C_\ell')
    plt.title("MWA, (%1.0e < flux < %1.0e) jy"%(mwa_jymin,mwa_jymax))
    
    plt.subplot(143)
    ispos = xspec>0
    isneg = xspec<0
    plt.loglog(lbins[ispos],xspec[ispos],'o')
    plt.loglog(lbins[isneg],-xspec[isneg],'v')
    plt.loglog(lbins,np.sqrt(pspec1*pspec2/bin_counts))
    plt.xlim([np.min(lbins),np.max(lbins)])
    plt.ylim([1.e-15,1.e-8])
    plt.xlabel('\ell')
    plt.ylabel('C_\ell')
    plt.title('cross spectrum')
    
    plt.subplot(144)
    xcor = xspec/np.sqrt(pspec1*pspec2)
    plt.errorbar(lbins,xcor,np.sqrt(.5*(1+xcor)/bin_counts))
    plt.ylim([-.1,.4])
    plt.xlim([0,np.max(lbins)])
    plt.xlabel('\ell')
    plt.ylabel('c_\ell')
    plt.title('cross correlation')
    plt.plot([0,5000],[0,0],'k-')
    
    plt.tight_layout()
    
 
    
def logloghist(plt,dat,bin_min,bin_max,num_bins,col):
    counts,bins = np.histogram(dat, 10.**np.linspace(np.log10(bin_min),np.log10(bin_max),num_bins+1))
    print(len(counts),len(bins))
    
    plt.loglog([bins[0],bins[0]],[.000001,counts[0]],col)
    for i in range(num_bins-1):
        plt.loglog([bins[i],bins[i+1],bins[i+1]],[counts[i],counts[i],counts[i+1]],col,lw=2)
    plt.loglog([bins[num_bins-1],bins[num_bins],bins[num_bins]],[counts[num_bins-1],counts[num_bins-1],.01],col,lw=2)
    
        
    plt.ylim([10,1.6*np.max(counts)])
    