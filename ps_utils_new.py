import numpy as np
from numpy import pi
import scipy.signal
from astropy.io import fits
#from scipy.optimize import curve_fit
import scipy.optimize

class IRImage:
    def __init__(self,fits_path,poly_order=4,run_poly_fit=False):
        print('loading '+fits_path)
        
        hdulist = fits.open(fits_path)
        self.header = hdulist[0].header
        self.full_ADU = hdulist[0].data
        self.n_full = self.full_ADU.shape[0]
        hdulist.close()
        
        self.MAGZPT = self.header['MAGZPT']
        self.exp_time_sec = self.header['EXPTIME']
        #self.filter_name = self.header['FILTER']
        self.dtheta_deg = float(np.abs(self.header['CD1_1']))
        self.dtheta_amin = 60.*self.dtheta_deg
        self.dtheta_rad = self.dtheta_deg*np.pi/180.
        self.domega_sr = self.dtheta_rad**2
        self.fov_deg = 5.
        self.n = int(np.round(self.fov_deg/self.dtheta_deg))
        
        self.ADU_to_kjy_per_sr = 3.631*(10**(-self.MAGZPT/2.5))/self.exp_time_sec/self.domega_sr
        self.full_kjy_per_sr = self.full_ADU * self.ADU_to_kjy_per_sr
        
        self.frame_ADU = self.crop_around_frame(self.full_ADU)
        self.ADU_to_rawADU = (1.86/60**2*np.pi/180.)**2/(self.domega_sr)
        self.frame_rawADU = self.frame_ADU*self.ADU_to_rawADU
        self.frame_kjy_per_sr = self.frame_ADU * self.ADU_to_kjy_per_sr
        
        if run_poly_fit: self.fit_2D_poly(poly_order)
        
    def crop_around_frame(self,full_img):
        x,y = np.where(full_img != 0)
        self.x_med, self.y_med = np.median(x),np.median(y)
        
        x1 = int(self.x_med-self.n/2+1)
        x2 = x1+self.n
        y1 = int(self.y_med-self.n/2+1)
        y2 = y1+self.n
    
        return full_img[x1:x2,y1:y2]
    
    def fit_2D_poly(self,order):
        def func(ind, *c):
            xi = 1.*(ind % self.n)/self.n-.5
            yi = 1.*(np.int32(ind/self.n))/self.n-.5

            xcoefs = c[1:1+order]
            ycoefs = c[1+order:]

            xpoly = 1 + np.sum([xcoefs[i]*(xi**(i+1)) for i in range(order)],axis=0)
            ypoly = 1 + np.sum([ycoefs[i]*(yi**(i+1)) for i in range(order)],axis=0)

            return xpoly*ypoly * c[0]  

        inds = np.arange(self.n**2)
        #popt, pcov = curve_fit(func, inds, self.frame_kjy_per_sr.flatten(), p0 = np.append(self.frame_kjy_per_sr.mean(),np.zeros(2*order)),method='trf')
        #np.sqrt(np.diag(pcov))

        chisq = lambda p0: np.sum((func(inds,*p0)-self.frame_rawADU.flatten())**2)
        bg = self.frame_rawADU.mean()
        print('bg = '+str(bg))
        poly_term_param_bounds = [(-.1,.1),(-.1,.1),(-.3,.3),(-.5,.5),(-.6,.6),(-1,1),(-2,2),(-2,2)]
        param_bounds = [(bg-10,bg+10)]+poly_term_param_bounds[:order] + poly_term_param_bounds[:order]
        print(param_bounds)

        result = scipy.optimize.differential_evolution(chisq, param_bounds, polish=True)
        popt = result['x']
        print(popt)
                               
        self.model_frame_rawADU = np.reshape(func(inds,*popt),(self.n,self.n))
        self.model_frame_kjy_per_sr = self.model_frame_rawADU / self.ADU_to_rawADU * self.ADU_to_kjy_per_sr
        self.res_frame_kjy_per_sr = self.frame_kjy_per_sr - self.model_frame_kjy_per_sr
        

def image2PS(image,nbins,lmax,backsub=False,hann=False,exclude_lx_or_ly_zero=False,use_res_image=False,use_model_image=False):
    n = image.n
    dang = image.dtheta_rad
    
    if use_res_image:
        img = image.res_frame_kjy_per_sr
    elif use_model_image:
        img = image.model_frame_kjy_per_sr
    else:
        img = image.frame_kjy_per_sr
    
    lvals = np.fft.fftfreq(n)*2*pi/dang
    dl = np.abs(lvals[1]-lvals[0])
    lx,ly = np.meshgrid(lvals,lvals)
    lmag  = np.sqrt(lx**2+ly**2)
    
    w2 = np.ones((n,n))
    norm = 1
    if hann:
        w = scipy.signal.hann(n)
        wx,wy = np.meshgrid(w,w)
        w2 = wx*wy
    norm = np.sqrt(np.mean(w2**2))
    
    sub = 0
    if backsub: sub = img.mean()
    print 'img.mean() = '+str(img.mean())
    img_ft = np.fft.fft2((img-sub)*w2)/norm
    
    lbinedges = np.linspace(0,lmax,nbins+1)
    lbincenters = .5*(lbinedges[0:nbins]+lbinedges[1:nbins+1])
    img_ft_binned = np.zeros(nbins)
    counts = np.zeros(nbins)
    for bini in range(nbins):
        inbin = (lmag>lbinedges[bini])&(lmag<lbinedges[bini+1])
        if exclude_lx_or_ly_zero: inbin &= (np.abs(lx)>3*dl) & (np.abs(ly)>3*dl)
        img_ft_binned[bini] = np.sum(np.abs(img_ft[inbin])**2)/np.sum(inbin)
        counts[bini] = np.sum(inbin)
    
    return lbincenters,img_ft_binned*(dang**2)/(n**2),counts

def make_hann(n):
    w = scipy.signal.hann(n)
    wx,wy = np.meshgrid(w,w)
    w2 = wx*wy
    return w2, np.sqrt(np.mean(w2**2))

def ir_and_radio_xspec(ir_image,ir_label,mwa_image,mwa_label,nbins,lmax):
    
    assert ir_image.dtheta_rad == mwa_image.dtheta_rad
    assert ir_image.n == mwa_image.n
   
    if ir_label == 'res':
        ir_img = ir_image.res_frame_kjy_per_sr
    else:
        ir_img = ir_image.frame_kjy_per_sr
        
    mwa_weights_img = mwa_image.weights_xx0
    if mwa_label == 'res':
        mwa_dirty_img = mwa_image.dirty_xx_u0-mwa_image.model_xx_u0
    else:
        mwa_dirty_img = mwa_image.dirty_xx_u0
    
    n,dang = ir_image.n, ir_image.dtheta_rad
    hann2D,hann2Drms = make_hann(n)
    lvals = np.fft.fftfreq(n)*2*pi/dang
    lx,ly = np.meshgrid(lvals,lvals)
    lmag  = np.sqrt(lx**2+ly**2)

    # FFT the (MWA dirty image) and (IR image)
    ir_ft = np.fft.fft2((ir_img-ir_img.mean())*hann2D)/hann2Drms
    mwa_dirty_ft = np.fft.fft2((mwa_dirty_img-mwa_dirty_img.mean()))
    mwa_weights_ft = np.abs(np.fft.fft2((mwa_weights_img-mwa_weights_img.mean())))
    
    lbinedges = np.linspace(0,lmax,nbins+1)
    lbincenters = .5*(lbinedges[0:nbins]+lbinedges[1:nbins+1])
    
    xspec_binned = np.zeros(nbins)
    mwaspec_binned = np.zeros(nbins)
    irspec_binned = np.zeros(nbins)
    
    bin_counts = np.zeros(nbins)
    bin_sum_weights = np.zeros(nbins)
    bin_sum_squared_weights = np.zeros(nbins)
    for bini in range(nbins):
        inbin = (lmag>lbinedges[bini])&(lmag<lbinedges[bini+1])
        bin_counts[bini] = np.sum(inbin)

        xspec_binned[bini] = np.sum(ir_ft[inbin]*np.conj(mwa_dirty_ft[inbin])*mwa_weights_ft[inbin])/np.sum(mwa_weights_ft[inbin])
        mwaspec_binned[bini] = np.sum(np.abs(mwa_dirty_ft[inbin])**2*mwa_weights_ft[inbin]**2)/np.sum(mwa_weights_ft[inbin]**2)
        irspec_binned[bini] = np.mean(np.abs(ir_ft[inbin])**2)
        
        bin_sum_weights[bini] = np.sum(mwa_weights_ft[inbin])
        bin_sum_squared_weights[bini] = np.sum(mwa_weights_ft[inbin]**2)
        
    pspec_norm = (dang**2)/(n**2)
    return lbincenters,irspec_binned*pspec_norm,mwaspec_binned*pspec_norm,xspec_binned*pspec_norm,bin_counts,bin_sum_weights,bin_sum_squared_weights
    
    
    
    