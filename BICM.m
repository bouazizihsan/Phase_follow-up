clc;
clear all;
warning off;
snrdb = [0:2:15];
snrLin = 10.^(snrdb/20);
sigma_p=0.128;
M=4;
mm=log2(M);
N = 50;
Ns = 500;

hMod = comm.RectangularQAMModulator ( 'ModulationOrder' , M , 'BitInput' , true,'NormalizationMethod','Average power');
hDemod = comm.RectangularQAMDemodulator ( 'ModulationOrder' , M,'BitOutput',1,'NormalizationMethod','Average power', 'DecisionMethod','Log-likelihood ratio');
pnoise = comm.PhaseNoise('Level',-50,'FrequencyOffset',20);

xCons=zeros(1,M);
 xCons = qammod(0:M-1,M,0, 'gray');
 power=sum(abs(xCons).^2)./length(xCons) ;
xCons= xCons./sqrt(power);


sign= randi([0,M-1],Ns,1);

Bin= de2bi(sign,'left-msb');
Bin=Bin';
Bin=Bin(:);
x= xCons(sign+1);
x=x.';
BinestBen=zeros(1,length(Bin));
llr=zeros(Ns,mm);
llr_ph=zeros(Ns,mm);
be=zeros(Ns,mm);

dataMod = step (hMod, Bin);


%  nf = modnorm(xCons,'peakpow',1);
%  x = nf*x;
%  max(x.*conj(x)) %pour verifier que l'energie =1
for s= 1:length(snrdb)
    s
    snr= snrLin(s);
    cont = 0;
    contM = 0;
    contBen = 0;
    count_ph=0;
    for h=1:N

        z  = complex(randn(size(x)),randn(size(x)))./ sqrt(2*snr);
        y  = x + z;
        yMod = dataMod + z;
        y_ph=exp(1j*sqrt(sigma_p)*randn(size(x))).*x +z ;
        
        dataOutM = step (hDemod, yMod);
        BinoutM = 1*(dataOutM<0) + 0*(dataOutM>=0);
        
        lBen = LLR_MQAM(M,y,xCons,snr);
        lBen = lBen';
        lBen = lBen(:);
        BinestBen = (lBen>0);
        llr_ph = LLR_MQAMphase(M,y_ph,xCons,snr,sigma_p);
        llr_ph = llr_ph';
        llr_ph = llr_ph(:);
        BinestBen_ph =  (llr_ph>0) ;
        contM = contM +  (length(Bin) - sum((Bin==BinoutM)));
        contBen = contBen +  (length(Bin) - sum((Bin==BinestBen)));
        count_ph=count_ph+(length(Bin)-sum(BinestBen_ph==Bin));
    end
    BErr_ph(s)=count_ph./(N*length(Bin))
    BErM(s) = contM/(length(Bin)*N);
    BErBEn(s) = contBen/(length(Bin)*N);
end

save rel5.mat




