clear all;
clc;

EbN0 = 5; %dB 
M = 16; % 2, 4, 16, 64 .... - QAM 
fading = 0; % fading on/off

% get M-QAM constellation and mapping 
h_mw = modem.qammod('M',M,'InputType', 'Bit','SymbolOrder','Gray'); 
constellation = h_mw.Constellation; 
mapping = h_mw.SymbolMapping;

% get variance of signal constellation (needed for the calculation of the 
% noise variance) 

sigma_a2 = sum(  abs( ( constellation - mean(constellation) ) ).^2/length(constellation)  );

% calculation of the noise variance 
EsN0 = 10*log10(10.^(EbN0/10) * log2(M)); 
sigma_n2 = sigma_a2./(10.^(EsN0./10));

llr_demod_obj = mex_llr_demod(constellation,mapping,'exact'); % generate object

for i = 1:100 % demodulate 100 frames of 1000 M-QAM Symbols 

% generate bitstream 
u = randi([0 1],1000*log2(M),1); 

% do modulation with matlab object 
x = modulate(h_mw,u); 

% generate fading factors 
if fading 
H = 1/sqrt(2)*randn(size(x)) + 1i * 1/sqrt(2) * randn(size(x)); 
else 
H = 1; 
end 

% calculate received signal 
r = x.*H + sqrt(0.5*sigma_n2)*randn(size(x)) + 1j * sqrt(0.5*sigma_n2)*randn(size(x)); 

% do Zero Forcing Equalization 
r_hat = r./H; 

% do demodulation (with variable NOISE_VAR sigma_n2./abs(H).^2) 
llr = calc_llr(llr_demod_obj, r_hat, sigma_n2./abs(H).^2); 
end