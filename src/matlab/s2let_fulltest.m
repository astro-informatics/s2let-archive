% flag_fulltest - Run all tests
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

% Main parameters
L = 32
N = 32
R = 10.0    

% Generate random 3D FLAG decomposition
flmn = zeros(N, L^2);
flmn = rand(size(flmn)) + sqrt(-1)*rand(size(flmn));
flmn = 2.*(flmn - (1+sqrt(-1))./2);

% Test complex synthesis - analysis exactness
f = flag_synthesis(flmn, L, N);
flmn_rec = flag_analysis(f, L, N);
flag_transform_error = max(max(abs(flmn-flmn_rec)))

% Impose reality on flms.
for en = 1:N
   for el = 0:L-1
      ind = el*el + el + 1;
      flmn(en,ind) = real(flmn(en,ind));
      for m = 1:el
         ind_pm = el*el + el + m + 1;
         ind_nm = el*el + el - m + 1;
         flmn(en,ind_nm) = (-1)^m * conj(flmn(en,ind_pm));
      end  
   end
end
% Test real synthesis - analysis exactness
f = flag_synthesis(flmn, L, N, 'reality', true);
flmn_rec = flag_analysis(f, L, N, 'reality', true);
flag_real_transform_error = max(max(abs(flmn-flmn_rec)))

% Generate random 1D SLAG decomposition
fn = rand(1,N);

% Test synthesis - analysis exactness
[f, nodes] = slag_synthesis(fn, N, 'R', R);
fn_rec = slag_analysis(f, N, R);
slag_transform_error = max(max(abs(fn-fn_rec)))
nodes2 = slag_sampling(N, R);
if (max(abs(nodes-nodes2))) ~= 0
    print('Problem with sampling scheme');
end
[f2, nodes] = slag_synthesis(fn, N, 'Nodes', nodes);
if max(max(abs(f-f2))) ~= 0
    print('Problem with transform when nodes is specified');
end

