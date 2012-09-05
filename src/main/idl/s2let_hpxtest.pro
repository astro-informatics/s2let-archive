pro s2let_hpxtest, B=B, L=L, J_min=J_min

if not keyword_set(B) then B = 7
if not keyword_set(L) then L = 192
if not keyword_set(J_min) then J_min = 2

loc = GETENV('S2LET')

if loc eq '' then begin

   print,'You must define environment variance S2LET (base directory of s2let)'

endif else begin

   file = loc + '/data/some_cmb_simu.fits'
   read_fits_map, file, f

   f_wav = s2let_axisym_hpx_wav_analysis(f, B, L, J_min)
   f_rec = s2let_axisym_hpx_wav_synthesis(f_wav, B, L, J_min)

   J = s2let_j_max(L, b)
   mollview, f_rec, title='Band-limited map'
   mollview, f_wav(*,0), title='Scaling map'
   for j=0, J-J_min do begin
      mollview, f_wav(*,j+1), title='Wavelet map '+strtrim(j,2)+' on '+strtrim(J-J_min+1,2)
   endfor

endelse

end
