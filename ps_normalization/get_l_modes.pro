function get_l_modes, d1, d2, pixx, pixy, $
                      pix=pix, $
                      verbose=verbose, $
                      ret_lx=ret_lx,$
                      ret_ly=ret_ly


;+
;
; Returns l-modes for dimensions [d1, d2]. Pixscale in arcmin/pix.
;
;-


if ~keyword_set(verbose) then verbose = 0


;;- Get fourier modes for d1
n = float(d1)
l = fltarr(n)
l1 = findgen(ceil(n/2))
l2 = (-1) * reverse( findgen(floor(n/2)) + 1)
l[0:(n_elements(l2)-1)] = l2
l[n_elements(l2):(n_elements(l)-1)] = l1
lx = l / (n / 2)

;;- Check l was filled properly.
if n_elements(l1) + n_elements(l2) ne n_elements(lx) then $
   message,'numel(l1) + numel(l2) != numel(lx)'


;;- Get fourier modes for d2
n = float(d2)
l = fltarr(n)
l1 = findgen(ceil(n/2))
l2 = (-1) * reverse( findgen(floor(n/2)) + 1)
l[0:(n_elements(l2)-1)] = l2
l[n_elements(l2):(n_elements(l)-1)] = l1
ly = l / (n / 2)
if n_elements(l1) + n_elements(l2) ne n_elements(ly) then $
   message,'numel(l1) + numel(l2) != numel(ly)'


;;- Make grid from l-modes
lx = shift(lx, ceil(d1/2.))
ly = shift(ly, ceil(d2/2.))
meshgrid, lx, ly, grid_lx, grid_ly


;;- l => 2pi/radian or /pix
IF keyword_set(pix) THEN BEGIN
   l = sqrt(grid_lx^2 + grid_ly^2) 
   message, 'ell (/pix)', /INF
ENDIF ELSE BEGIN
   grid_lx *= 2. * !pi / (pixx / 3600. * !pi / 180.) / 2.
   grid_ly *= 2. * !pi / (pixy / 3600. * !pi / 180.) / 2.
  ; extra factor of 1/2 for nyquist limit, bro.
                                ;grid_lx /= (pixscale * !pi / 180.)
                                ;grid_ly /= (pixscale * !pi / 180.)
   l = sqrt(grid_lx^2 + grid_ly^2)
   if verbose then message, 'ell [2 pi/theta] in units (1/rad)', /INF
ENDELSE

if keyword_set(ret_lx) then return, lx
if keyword_set(ret_ly) then return, ly

RETURN, l

END
