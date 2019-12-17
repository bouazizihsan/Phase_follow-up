% 4 QAM
clear all;
clc;
snrdb =[0:2:10];
snrlin= 10.^(snrdb./10);
M=16;
m=log2(M);
Npaquet=100;
sigma_p =0.128;

sig= randi([0,1],1,1000*m);

sign=reshape(sig,[m,length(sig)./m])';

Sym=bi2de(sign, 'left-msb');
x= qammod(Sym,M,0,'gray');

Binest=zeros(1,length(sig));
Binest_ph=zeros(1,length(sig));
xCons=zeros(1,M);
for id=1:M
    xCons(id)=qammod(id-1,M,0,'gray'); %save constellation
end

Ns = length(x);
llr=zeros(Ns,m);
llr=zeros(Ns,m);
llr=zeros(Ns,m);
% refconst = qammod(0:M-1,M);
% res=refconst - xCons;

% nf = modnorm(xCons,'peakpow',1);
% sigMod = nf*sigMod;
% max(sigMod.*conj(sigMod))

for i=1:length(snrdb)
    i
    snr= snrlin(i);
    count_ph=0;
    count=0;
   count_ph_BLT=0;
 for s=1:Npaquet 
     
    z=complex(rand(size(x)),rand(size(x)))./sqrt(2.*snr);
    y_ph=exp(1j*randn(size(x))*sigma_p).*x +z ;
     y=x +z ;
     
      llr = LLR_MQAM(M,y,xCons,snr);
       llr = llr';
        llr = llr(:);
        Binest =(llr>0) ;
        
        [llr_ph_LT,llr_ph_BLT ]= LLR_MQAMphase(M,y_ph,xCons,snr,sigma_p);
        llr_ph_LT = llr_ph_LT';
        llr_ph_LT = llr_ph_LT(:);
        Binest_ph =  (llr_ph_LT>0) ;
        
        llr_ph_BLT = llr_ph_BLT';
        llr_ph_BLT = llr_ph_BLT(:);
        Binest_ph_BLT =  (llr_ph_BLT>0) ;
        
   count   =count+(length(sig)-sum(Binest==sig'));
   
   count_ph=count_ph+(length(sig)-sum(Binest_ph==sig'));
   
   count_ph_BLT=count_ph_BLT+(length(sig)-sum(Binest_ph_BLT==sig'));
   
end

BErr(i)=count./(Npaquet.*length(sig));
BErr_ph(i)=count_ph./(Npaquet.*length(sig));
BErr_ph_BLT(i)=count_ph_BLT./(Npaquet.*length(sig));

end

save phase.mat
