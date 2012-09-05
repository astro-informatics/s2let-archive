function s2let_axisym_hpx_wav_analysis, f, B, L, J_min
; S2LET package
; Copyright (C) 2012 
; Boris Leistedt & Jason McEwen

loc = GETENV('S2LET')

if loc eq '' then begin

   print,'You must define environment variance S2LET (base directory of s2let)'
   return, 0

endif else begin

   soname = loc + '/lib/c/libs2let.so'

   nside = round(sqrt( (size(f))(1) / 12.0 ))
   npix = nside2npix(nside)

   J = s2let_j_max(L, B)
   help,nside,L, B, J_min, npix, J

   f_wav_vec = dblarr((J+1-J_min)*npix)
   f_scal = dblarr(npix)

   r = call_external(soname, 's2let_idl_axisym_hpx_wav_analysis_real', f_wav_vec, f_scal, f, nside, B, L, J_min, /CDECL)
   
   f_wav = dblarr(npix, J+2-J_min)
   f_wav(0:npix-1, 0) = f_scal(0:npix-1)
   for j=0, J-J_min do begin
      f_wav(0:npix-1, j+1) = f_wav_vec( j*npix : (j+1)*npix-1 )
   endfor

   return, f_wav

endelse

end
