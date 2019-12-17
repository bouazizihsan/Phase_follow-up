clear all;
clc;

EbN0 = 5; %dB 
M = 16; % 2, 4, 16, 64 .... - QAM 
fading = 1; % fading on/off
N=100;
% get M-QAM constellation and mapping 
h_mw = modem.qammod('M',M,'InputType', 'Bit','SymbolOrder','Gray'); 
constellation = h_mw.Constellation; 
mapping = h_mw.SymbolMapping;

% get variance of signal constellation (needed for the calculation of the 
% noise variance) 

% sigma_a2 = sum(  abs( ( constellation - mean(constellation) ) ).^2/length(constellation)  );

% calculation of the noise variance 
% EsN0 = 10*log10(10.^(EbN0/10) * log2(M)); 
% sigma_n2 = sigma_a2./(10.^(EsN0./10));

llr_demod_obj = mex_llr_demod(constellation,mapping,'exact'); % generate object

snrdb = [0:2:15];
snrLin = 10.^(snrdb/20);

% generate bitstream 
u = randi([0 1],1000*log2(M),1); 

% do modulation with matlab object 
x = modulate(h_mw,u); 
if fading 
H = 1/sqrt(2)*randn(size(x)) + 1i * 1/sqrt(2) * randn(size(x)); 
else 
H = 1; 
end 

for s= 1:length(snrdb)
    s
    snr= snrLin(s);
    cont = 0;
for i = 1:N 

% calculate received signal 
r = x.*H + sqrt(0.5*snr)*randn(size(x)) + 1j * sqrt(0.5*snr)*randn(size(x)); 

% do Zero Forcing Equalization 
r_hat = r./H; 

% do demodulation (with variable NOISE_VAR sigma_n2./abs(H).^2) 
llr = calc_llr(llr_demod_obj, r_hat, snr./abs(H).^2); 

llr = llr';
        llr = llr(:);
        BinestBen =  1*(llr>0) + 0*(llr<=0);

        cont = cont +  (length(u) - sum((u==BinestBen)));
end

 BErBEn(s) = cont/(length(u)*N)

end