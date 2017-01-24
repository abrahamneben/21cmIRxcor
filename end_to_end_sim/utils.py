import numpy as np
import scipy.signal

def make_hann_3D(n):
    w = scipy.signal.hann(n)
    wx,wy,wz = np.meshgrid(w,w,w)
    w3 = wx*wy*wz
    return w3, np.sqrt(np.mean(w3**2))

def make_hann_2D(n):
    w = scipy.signal.hann(n)
    wx,wy = np.meshgrid(w,w)
    w2 = wx*wy
    return w2, np.sqrt(np.mean(w2**2))

def calc_3D_xspec(cube1,cube2,dx,nbins,kmin,kmax,uselogbins=False,usehann=True):
   
    N = cube1.shape[0]
    V = (1.*N*dx)**3
    
    hann3D,hann3Drms = np.ones((N,N,N)),1
    if usehann:
        hann3D,hann3Drms = make_hann_3D(N)
    
    dk = 2.*np.pi/(N*dx)
    kvals = np.fft.fftfreq(N,1)*N*dk
    kxgrid,kygrid,kzgrid = np.meshgrid(kvals,kvals,kvals)
    kmaggrid = np.sqrt(kxgrid**2+kygrid**2+kzgrid**2)
    print(np.max(kmaggrid))

    cube1_ft = dx**3*np.fft.fftn((cube1-cube1.mean())*hann3D)/hann3Drms
    cube2_ft = dx**3*np.fft.fftn((cube2-cube2.mean())*hann3D)/hann3Drms
    
    if uselogbins: 
        kbinedges = 10.**np.linspace(np.log10(kmin),np.log10(kmax),nbins+1)
    else:
        kbinedges = np.linspace(kmin,kmax,nbins+1)
    kbincenters = .5*(kbinedges[0:nbins]+kbinedges[1:nbins+1])
    
    xspec_binned = np.zeros(nbins)
    pspec1_binned = np.zeros(nbins)
    pspec2_binned = np.zeros(nbins)

    bin_counts = np.zeros(nbins)
    for bini in range(nbins):
        inbin = (kmaggrid>kbinedges[bini])&(kmaggrid<kbinedges[bini+1])
        bin_counts[bini] = np.sum(inbin)

        xspec_binned[bini] = np.real(np.mean(cube1_ft[inbin]*np.conj(cube2_ft[inbin])))
        pspec1_binned[bini] = np.mean(np.abs(cube1_ft[inbin])**2)
        pspec2_binned[bini] = np.mean(np.abs(cube2_ft[inbin])**2)
                
    return kbincenters,xspec_binned/V,pspec1_binned/V,pspec2_binned/V,bin_counts
    
def calc_2D_xspec(img1,img2,dtheta_rad,nbins,lmin,lmax,uselogbins=False,usehann=True):

    N = img1.shape[0]
    hann2D,hann2Drms = make_hann_2D(N)
    lvals = np.fft.fftfreq(N)*2*np.pi/dtheta_rad
    dl = np.abs(lvals[1]-lvals[0])
    lx,ly = np.meshgrid(lvals,lvals)
    lmag  = np.sqrt(lx**2+ly**2)
    print('dl = %d\nlmax = %d\n'%(dl,np.max(lmag)))

    # FFT the (MWA dirty image) and (IR image)
    img1_ft = np.fft.fft2((img1-img1.mean())*hann2D)/hann2Drms
    img2_ft = np.fft.fft2((img2-img2.mean())*hann2D)/hann2Drms
    
    if uselogbins: 
        lbinedges = 10.**np.linspace(np.log10(lmin),np.log10(lmax),nbins+1)
    else:
        lbinedges = np.linspace(lmin,lmax,nbins+1)
    lbincenters = .5*(lbinedges[0:nbins]+lbinedges[1:nbins+1])
    
    xspec_binned = np.zeros(nbins)
    pspec1_binned = np.zeros(nbins)
    pspec2_binned = np.zeros(nbins)

    bin_counts = np.zeros(nbins)
    for bini in range(nbins):
        inbin = (lmag>lbinedges[bini])&(lmag<lbinedges[bini+1])
        bin_counts[bini] = np.sum(inbin)

        xspec_binned[bini] = np.mean(img1_ft[inbin]*np.conj(img2_ft[inbin]))
        pspec1_binned[bini] = np.mean(np.abs(img1_ft[inbin])**2)
        pspec2_binned[bini] = np.mean(np.abs(img2_ft[inbin])**2)

    pspec_norm = (dtheta_rad**2)/(N**2)
    return lbincenters,xspec_binned*pspec_norm,pspec1_binned*pspec_norm,pspec2_binned*pspec_norm,bin_counts
    
    