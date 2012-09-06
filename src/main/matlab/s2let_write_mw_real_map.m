function s2let_write_mw_real_map(f, file)

sz = size(f);
L = min([sz(1), sz(2)]);

mwmap = [];
for l = 1:L
    mwmap = [mwmap f(l,:)];
end

fitswrite(mwmap, file, 'L', L);

end