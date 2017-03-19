FUNCTION get_raw_spec, array, array2, pixx, pixy,$
                       fname=fname,$                    
                       nbins=nbins, $
                       linear=linear, $                       
                       pad=pad,$
                       ibins=ibins,$
                       keep_all=keep_all,$
                       verbose=verbose,$
                       w=w, $
                       fft_=fft_
  

;+
;
; To compute the raw angular power spectrum of an image(s).
; pixsize should be in arcsec/pix.
;
; Ketron 12/2012
;
;-


COMPILE_OPT HIDDEN
tstart = systime(/sec)


dim1 = (size(array))[1]
dim2 = (size(array))[2]


if ~keyword_set(nbins) then nbins = 42 
if ~keyword_set(w) then w = fltarr(dim1, dim2) + 1.0
if n_elements(pixy) eq 0 then pixy = pixx
if keyword_set(ibins) then nbins = n_elements(ibins)
input_nb = nbins

;w[2400:2560,1510] = 0.0
;w[2475,1440:1561] = 0.0

;;- Second array size
if ( (size(array2))[1] ne dim1 ) or $
   ( (size(array2))[2] ne dim2 ) then $
      message, "Arrays must be the same size."

if keyword_set(pad) then begin
   ;;- Change the number of pixels slightly, for FFT speed.
;;- FFT running time is "roughly proportional to the 
;; total number of points in Array times the sum of 
;; its prime factors." -- edit array sizes.
   old_dim1 = dim1 & old_dim2 = dim2

   factor, dim1, prime1, count1, /quiet
   n_prime1 = n_elements(prime1)
   largest_prime1 = prime1[n_prime1 - 1]
   WHILE prime1[n_prime1-1] GT 15 DO BEGIN
      dim1 -= 1
      factor, dim1, prime1, count1, /quiet
      n_prime1 = n_elements(prime1)
   ENDWHILE

;;- Second dimension
   factor, dim2, prime2, count2, /quiet
   n_prime2 = n_elements(prime2)
   largest_prime2 = prime2[n_prime2 - 1]
   WHILE prime2[n_prime2-1] GT 15 DO BEGIN
      dim2 -= 1
      factor, dim2, prime2, count2, /quiet
      n_prime2 = n_elements(prime2)
   ENDWHILE
   
   d1 = dim1 - old_dim1
   d2 = dim2 - old_dim2

;;- Edit crop sizes.
   message, 'removing ' + strim(d1) + ' pixel(s) to the first dimension'+$
   'and ' + strim(d2) + ' pixel(s) to the second dimension.', /INF
   x1 = dim1 - d1   
   y1 = dim2 - d2  

   array = array[0:dim1-1,0:dim2-1]
   array2 = array2[0:dim1-1,0:dim2-1]

endif



;;- Check for NaNs -- this will screw up the FFT
nan = where(finite(array) eq 0)
if nan[0] NE -1 then message, 'NaNs in the array need to be ' + $
                              'incorporated into the mask.'



;;- Need to normalize (IDL returns 1/N*fft) and multiply by 
;; sampling interval in radians.
dim1 = (size(array))[1]
dim2 = (size(array))[2]
ft2 = w * real_part(fft(array)*conj(fft(array2))) * $
      (!pi*pixx/180./3600.) * (!pi*pixy/180./3600.) * (float(dim1)*float(dim2))

if keyword_set(fft_) then begin
   ft2 = shift(ft2, dim1/2, dim2/2)
   return, ft2
endif


;;- Binning
l = get_l_modes(dim1, dim2, pixx, pixy)
;message, 'Forcing max ell bin at 1e4',/inf 
if keyword_set(linear) then bins = linbins(nbins, 900., 6.55e6) else $
   bins = logbins(nbins, min(l[where(l gt 0)]), max(l));1e4
;stop
;bins = [1534.0909, 64800.0, 83113.942, 129515.23 ]
;nbins=3
;stop

;min(l[where(l gt 0)])-100.
if keyword_set(ibins) then bins = ibins


Cl = fltarr(nbins-1)
lbins = fltarr(nbins-1)
Nl = fltarr(nbins-1) ;total number of fourier elements per bin.
bin_edges = fltarr(nbins*2 + 1)
j = 0
for i = 0, n_elements(Cl) - 1 do begin
   inside = where((l ge bins[i]) and (l le bins[i+1]), Nbeans)
   Nl[i] = Nbeans
   if inside[0] ne -1 then begin
      ft2_in = ft2[inside]
      ft2_in = ft2_in[where(ft2_in ne 0)]
      w_in = w[inside]
      
;      if min(ft2_in) lt 0 then begin
;         stop
;         ft2_in[where(ft2_in eq min(ft2_in))] *= -1
;      endif

      cl[i] = total(ft2_in) / total(w_in)

      lbins[i] = mean([bins[i], bins[i+1]])
      bin_edges[j]= bins[i]
      bin_edges[j+1] = bins[i+1]
      j+=2
      ;print, 'Ninside = '+ strim(Nbeans)+ ' at '+strim(bins[i]) + ' and '+strim(bins[i+1])
   endif else begin
      ;print, 'Ninside = 0 at '+strim(bins[i]) + ' and '+strim(bins[i+1])
   endelse
endfor
nonz = where((cl ne 0)); and (Nl ge 10)) ;<---- REMOVING LOW ELL!!!!!!
cl = cl[nonz]
l = lbins[nonz]
Nl = Nl[nonz]
if (min(nonz))[0] ne 0 then bin_edges[0:min(nonz)*2-1]=0
bin_edges = bin_edges[where(bin_edges ne 0)]

cl_struct = { l:fltarr(n_elements(l)),$
              Cl:fltarr(n_elements(Cl)),$
              bin_edges:fltarr(n_elements(bin_edges)),$
              dcl:fltarr(n_elements(cl)),$
              Nl:fltarr(n_elements(l))}

cl_struct.l = l
cl_struct.Cl = Cl
cl_struct.bin_edges = bin_edges
cl_struct.Nl = Nl


;;- Cosmic variance -> dCl
fishdoggin = where((array ne 0) and (array2 ne 0), Ngood)
pix = sqrt(pixx^2 + pixy^2)
fsky = double(Ngood) * (pix/3600.)^2 * (!pi/180.)^2 ;sr
delta_l = fltarr(n_elements(l))
j = 0
for i = 0, n_elements(bin_edges) - 2, 2 do begin
   delta_l[j] = bin_edges[i+1]-bin_edges[i]
   j+=1
endfor
fsky /= (4*!pi) ;fraction of total sky. 
cosmic_variance = sqrt( 2 / (((2 * l) + 1) * fsky * delta_l))
cl_struct.dcl = cosmic_variance ;* cl



t_tot = (systime(/sec) - tstart) / 60.
if keyword_set(verbose) then message, strim(t_tot, l=4) + ' minutes.', /INF 
nbins = input_nb
return, cl_struct


end
