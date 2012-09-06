function [healpixmap, nside]= s2let_read_hpx_real_map(file)

datacell = fitsread(file,'binarytable');
data = datacell{1};
sz = size(data);

healpixmap = [];
for col = 1:sz(1)
    healpixmap = [healpixmap data(col,:)];
end

nside = sqrt((sz(1) * sz(2)) / 12);

end