import numpy as np
from numpy import pi,linspace,array,fft,abs,mean,floor,zeros,sum,sqrt,sort,savez
import healpy as hp
from scipy.io import readsav

c=299792458.
kb=1.3806488e-23
jy=1.e-26

# two use cases:
# (1) load a 10deg FOV centered at EOR0, apply uniform weighting to the whole field
# (2) load a 15deg FOV centered at EOR0, apply uniform weighting to a sub field of size n_crop, centered at x_med,y_med
class MWAImage:
    def __init__(self,freq_averaged_cubedat0,freq_averaged_cubedat1,dtheta_amin,n,crop_before_uniform_weighting_params=None):
        print('initializing MWAImage object')
        
        self.theta0_rad = -27.*pi/180.
        self.phi0_rad = 0
        
        self.dtheta_amin = 1.*dtheta_amin
        self.dtheta_rad = self.dtheta_amin/60*np.pi/180.
        self.n = n
        self.angvals_rad = np.linspace(-self.dtheta_rad*self.n/2,self.dtheta_rad*((self.n-1.)/2.),self.n)        
        self.freqs = np.linspace(168.e6,199.e6,384)
        
        self.dirty_xx0,self.counts_cube,self.weights_xx0,self.model_xx0,self.beam_squared_xx0 = grid_fhd_cubes(\
                            freq_averaged_cubedat0, self.angvals_rad, self.phi0_rad, self.theta0_rad)
        self.dirty_xx1,self.counts_cube,self.weights_xx1,self.model_xx1,self.beam_squared_xx1 = grid_fhd_cubes(\
                            freq_averaged_cubedat1, self.angvals_rad, self.phi0_rad, self.theta0_rad)

        self.x_psf_cent,self.y_psf_cent = argmax2D(self.weights_xx0)
        print("psf is centered at (x,y) = (%d,%d)"%(self.x_psf_cent, self.y_psf_cent))
        
        if crop_before_uniform_weighting_params != None:
            print('cropping before applying uniform weighting')
            self.n, self.x_med, self.y_med = crop_before_uniform_weighting_params
            
            self.dirty_xx0 = crop(self.dirty_xx0, self.n, self.x_med, self.y_med)
            self.dirty_xx1 = crop(self.dirty_xx1, self.n, self.x_med, self.y_med)
            self.model_xx0 = crop(self.model_xx0, self.n, self.x_med, self.y_med)
            self.model_xx1 = crop(self.model_xx1, self.n, self.x_med, self.y_med)
            
            self.weights_xx0 = crop(self.weights_xx0, self.n, self.x_psf_cent, self.y_psf_cent)
            self.weights_xx1 = crop(self.weights_xx1, self.n, self.x_psf_cent, self.y_psf_cent)
        

        self.dirty_xx_u0 = apply_uniform_weighting(self.dirty_xx0, self.weights_xx0, self.n,self.dtheta_rad,self.freqs)
        self.model_xx_u0 = apply_uniform_weighting(self.model_xx0, self.weights_xx0, self.n,self.dtheta_rad,self.freqs)

        self.dirty_xx_u1 = apply_uniform_weighting(self.dirty_xx1, self.weights_xx1, self.n,self.dtheta_rad,self.freqs)
        self.model_xx_u1 = apply_uniform_weighting(self.model_xx1, self.weights_xx1, self.n,self.dtheta_rad,self.freqs)

def load_freq_averaged_odd_even_cubedat(fhdcubesroot,fhdlabel,pol='XX'):
    return load_and_freq_average_cube(fhdcubesroot+fhdlabel+'_odd_cube'+pol+'.sav'),\
           load_and_freq_average_cube(fhdcubesroot+fhdlabel+'_even_cube'+pol+'.sav')
        
def load_and_freq_average_cube(cubepath):
    cubedat = readsav(cubepath)
    cubedat['dirty_cube'] = np.mean(cubedat['dirty_cube'],axis=0)
    cubedat['model_cube'] = np.mean(cubedat['model_cube'],axis=0)
    cubedat['variance_cube'] = np.mean(cubedat['variance_cube'],axis=0)
    cubedat['weights_cube'] = np.mean(cubedat['weights_cube'],axis=0)
    cubedat['res_cube'] = np.mean(cubedat['res_cube'],axis=0)
    cubedat['beam_squared_cube'] = np.mean(cubedat['beam_squared_cube'],axis=0)
    return cubedat
        
def argmax2D(img):
    x,y = np.where(img == np.max(img))
    return x[0],y[0]
        
def crop(img,n,x_med,y_med):
    print("cropping to (%d,%d)"%(x_med,y_med))
    x1 = int(x_med-n/2+1)
    x2 = x1+n
    y1 = int(y_med-n/2+1)
    y2 = y1+n

    return img[x1:x2,y1:y2]
        
def healpix_cube_to_xy_cube(healpixi_cube,hpix_inds,nside,angvals,phi0,theta0,return_counts_cube=False):
    print('gridding healpix')
    minangval = np.min(angvals)
    dangval = angvals[1]-angvals[0]
    n = len(angvals)

    v = hp.pix2vec(nside,hpix_inds)
    #vrot = hp.rotator.rotateVector(v,euler_matrix=euler_matrix_new(-phi0,0,0,/zyx))
    vrot = hp.rotator.rotateVector(rotmat=hp.rotator.euler_matrix_new(0,-(theta0-pi/2),0,ZYX=True),vec=v)

    x_indices = np.int32((vrot[0,:]-minangval)/dangval)
    y_indices = np.int32((vrot[1,:]-minangval)/dangval)

    xy_cube = zeros((n,n))
    counts_cube = zeros((n,n))
    for xi in range(n):
        for yi in range(n):
            matchingpixels = (x_indices==xi)*(y_indices==yi)
            xy_cube[xi,yi] = mean(healpixi_cube[matchingpixels])
            counts_cube[xi,yi] = sum(matchingpixels)
    if return_counts_cube == False: return xy_cube
    return xy_cube,counts_cube
    
def grid_fhd_cubes(cubedat,angvals,phi0,theta0):
    hpx_inds = cubedat['hpx_inds']
    nside = 1024

    def gridcube(hpix_cube,return_counts_cube=False):
        return healpix_cube_to_xy_cube(hpix_cube,hpx_inds,nside,angvals,phi0,theta0,return_counts_cube)

    dirty_xx, counts_cube = gridcube(cubedat['dirty_cube'],return_counts_cube=True)

    return dirty_xx,counts_cube,gridcube(cubedat['weights_cube']),\
           gridcube(cubedat['model_cube']),gridcube(cubedat['beam_squared_cube'])
        
    
def apply_uniform_weighting(cube_xy,weights_xy,n,dtheta_rad,freqs):
    print('applying uniform weighting')
    cube_uv = fft_xy2uv(cube_xy, n, dtheta_rad)
    weights_uv = abs(fft_xy2uv(weights_xy, n, dtheta_rad))

    cube_uv_u = cube_uv / weights_uv
    cube_uniform_xy = fft_uv2xy(cube_uv_u, n, dtheta_rad)
    cube_uniform_kelvin_xy = np.real(cube_uniform_xy) *jy*c**2/(2*kb*mean(freqs)**2)

    return cube_uniform_kelvin_xy


def fft_uv2xy(cube_uv, n,dtheta):
    du = 1./(n*dtheta)
    cube_xy = fft.fft2(cube_uv)*du**2
    return cube_xy
def fft_xy2uv(cube_xy, n,dtheta):
    cube_uv = fft.ifft2(cube_xy)*n**2*dtheta**2
    return cube_uv




def img2PS(img,weight_img,pixsize_rad,nbins,lmax,backsub=False,hann=False):
    n = img.shape[0]
    dang = pixsize_rad
    lvals = np.fft.fftfreq(n)*2*pi/dang
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
    #print 'img.mean() = '+str(img.mean())
    img_ft = np.fft.fft2((img-sub)*w2)/norm
    weight_ft = np.abs(np.fft.fft2(weight_img*w2))
    
    lbinedges = np.linspace(1,lmax,nbins+1)
    lbincenters = .5*(lbinedges[0:nbins]+lbinedges[1:nbins+1])
    img_ft_binned = np.zeros(nbins)
    img_ft_binned_weighted = np.zeros(nbins)
    bin_counts = np.zeros(nbins)
    bin_weights = np.zeros(nbins)
    for bini in range(nbins):
        inbin = np.logical_and(lmag>lbinedges[bini],lmag<lbinedges[bini+1])
        # img_ft_binned[bini] = np.sum(np.abs(img_ft[inbin])**2)/np.sum(inbin)
        bin_counts[bini] = np.sum(inbin)

        # weight_ft prop to observation time in that cell,
        # meaning it is inv prop to vavriance in that cell
        # so we should directly use it to weight the sum
        img_ft_binned_weighted[bini] = np.sum(np.abs(img_ft[inbin])**2*weight_ft[inbin]**2)/np.sum(weight_ft[inbin]**2)
        bin_weights[bini] = np.sum(weight_ft[inbin]**2)
    
    return lbincenters,img_ft_binned_weighted*(dang**2)/(n**2),bin_counts,bin_weights
