function llr = calc_llr(r,sigma_n2,constellation,mapping)

%CALC_LLR LLR-Value Computation
%   LLR = CALC_LLR(R, NOISE_VAR, CONSTELLATION, MAPPING) calculates the
%   LLR-values LLR based on the received signal R. The vectors CONSTELLATION 
%   and MAPPING are usually generated by using MATLAB Modem Modulation 
%   Modem.<Type> and contain the complex or real constellations points and 
%   the bit mapping. NOISE_VAR denotes the noise variance of the received
%   signal which must be available or at least an estimate for it.
%
%       Input:      R                       [Nx1]
%                   NOISE_VAR               [1x1]
%                   CONSTELLATION           [1xM]
%                   MAPPING                 [1xM]
%       Output:     LLR                     [Nx1]
%
%   Author: Bernhard Schmidt
%
% Copyright 2010 by Bernhard Schmidt
% Permission is granted to use this program/data for educational/research only

mod_idx = log2(length(constellation));

dist = zeros(2^mod_idx,length(r));

for idx = 1:(2^mod_idx) 
        dist(idx,:) = (real(r) - (real(constellation(idx)))).^2 + (imag(r) - (imag(constellation(idx)))).^2;
end

exp_mat = exp(-1./dist*sigma_n2);
llr = zeros(mod_idx,length(r));

for idx = 1:mod_idx
    llr(mod_idx-idx+1,:) = log( sum( exp_mat( (bitget( mapping,idx )==0),:) ,1) ) - log(sum(exp_mat((bitget(mapping,idx)==1),:),1));
end
llr=llr(:);