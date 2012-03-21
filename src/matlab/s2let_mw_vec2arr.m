function arr = s2let_mw_vec2arr(vec)

sz = size(vec);

if sz(1) == 1   % array 1 x L*(2L-1)
    delta = sqrt(1 + 8*(sz(2)));
    L = ( 1 + delta ) / 4;
else            % array L*(2L-1) x 1
    delta = sqrt(1 + 8*(sz(1)));
    L = ( 1 + delta ) / 4;
end

arr = zeros(L, 2*L-1);

for t = 1:L
    for p = 1:2*L-1
        arr(t,p) = vec((t-1)*(2*L-1)+p);
    end
end


end