function s2let_write_hpx_real_map(f, file)

sz = size(f);
szb = max([sz(1), sz(2)]);
nside = floor(sqrt(szb/12.0));

fitswrite(f, file, 'NSIDE', nside);

end