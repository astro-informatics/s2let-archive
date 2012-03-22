function error_on_axisym_tilling = s2let_check_axisym_tilling(kappa, kappa0, L, J)

% s2let_check_axisym_tilling - Checks exactness of the tilling.
% -- Axisymmetric wavelets on the sphere.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

identity = kappa0.^2;
for j=0:J
    identity(1,:) = identity(1,:) + kappa(j+1,:).^2;
end

error_on_axisym_tilling = 0;
for l=1:L
    error_on_axisym_tilling = error_on_axisym_tilling + identity(1,l) - 1.0; 
end

end