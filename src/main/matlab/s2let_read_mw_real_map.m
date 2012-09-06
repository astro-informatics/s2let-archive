function [mwmaparr, L]= s2let_read_mw_real_map(file)

datacell = fitsread(file,'binarytable');
data = datacell{1};
sz = size(data);

mwmap = [];
for col = 1:sz(1)
    mwmap = [mwmap data(col,:)];
end

mwmaparr = s2let_mw_vec2arr(mwmap);
sz = size(mwmaparr);
L = sz(1);

end