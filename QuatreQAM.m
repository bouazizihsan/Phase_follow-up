clc;
clear all;
warning off;
snrdb = [0:1:20];
snrLin = 10.^(snrdb/20);
Nb = 2^8;
M=4;
N = 100;


Sym=randi([0,3],1,Nb);
%
% x  = map4QAM(Bin);
% re = reshape(Bin,[2,length(Bin)/2])';
% be =bi2de(re, 'left-msb');
x = constellation(Sym);
xMatlab = qammod(Sym,M);

for s= 1:length(snrdb)
    s
    snr= snrLin(s);
    cont = 0;
    cont2 = 0;
    for h=1:N
        z  = complex(randn(size(x)),randn(size(x)))./ sqrt(2*snr);
        y  = x + z;
        yMatlab = xMatlab + z;
        [xDec,Binh] = dem(y,M);
        SymMatlab =  qamdemod(yMatlab,M);
        %      BinHMatlab = reshape(de2bi(matlabD,'left-msb')',[1,length(Bin)]);
        cont = cont +  (length(Sym) - sum((Sym==xDec)));
        cont2 = cont2 +  (length(Sym) - sum((Sym==SymMatlab)));
    end
    
    SyEr(s) = cont/(length(Sym)*N);
    SyErM(s) = cont2/(length(Sym)*N);
end

save rel5.mat




