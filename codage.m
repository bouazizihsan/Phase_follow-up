clc;
clear all;
warning off;
snrdb = [-5:1:20];
snrLin = 10.^(snrdb/20);
Nb = 2^3;
M=4;
N=7;
K=4;
Imax=2;
mm=log2(M);

 G = [1 0 0 0 1 0 1; 0 1 0 0 1 1 1;0 0 1 0 1 1 0; 0 0 0 1 0 1 1];
 H =[1 1 1 0 1 0 0; 0 1 1 1 0 1 0; 1 1 0 1 0 0 1];

%  G =[1 0 1; 0 1 1]; % C(3,2), size(G)= 2 3
%   H =[1 1 1];
% rate =4/7 .* 1/2 = 2/7
msg=randi([0 1],1,104);
re = reshape(msg,[K,length(msg)/K])';
msgcoded= re*G;
msgcoded=mod(msgcoded,2);
% ree=re';
% ree=ree(:);
Bink=msgcoded';
Binaire=Bink(:); 
Bin=Binaire';
msgcod=reshape(Bin,[2,length(Bin)/2])';
Sym=bi2de(msgcod, 'left-msb');
x = qammod(Sym,16,0,'gray');

Npaquet = 50;
% bb=msgcoded';
Ns = length(x);
be=zeros(Ns,mm);
% dataMod = step (hMod, msg');
xCons=zeros(1,M);
for id=1:M
    xCons(id)=qammod(id-1,M,0,'gray'); %save constellation
end
for s= 1:length(snrdb)
  s
    snr= snrLin(s);
    cont = 0;
    contt = 0;
    
    for h=1:Npaquet
        z  = complex(randn(size(x)),randn(size(x)))./ sqrt(2*snr);
        y  = x + z;
        
        llr = LLR_MQAM(M,y,xCons,snr);
        
        llr_bits=llr';
        llr_bits=llr_bits(:);
        llr_bits=llr_bits';
         llr=reshape(llr_bits,[N,length(llr_bits)./N])';
         for comp=1:size(llr,1)
              xMP(comp,:)=MsgPassAlg(llr(comp,:),H,Imax);
         end
        llrmp=xMP';
        llrmp= llrmp(:);
        Bin =  1*(llrmp>0) + 0*(llrmp<=0);
        inter=Bin';
        Bines=reshape(inter,[N,length(inter)/N])';
%         Syndrom=Bines*H'; 
%          Syndrom=mod(Syndrom,2);
        bindec = Bines';
       bindec = bindec(:);
       
        llrT = llr';
        llrT = llrT(:);
        BinestBen =  1*(llrT>0) + 0*(llrT<=0);
        aux=BinestBen';
        Binest=reshape(aux,[N,length(aux)/N])';
%         Synd=Binest*H'; 
%         Synd=mod(Synd,2);
       bindecoded= Binest; 
       bindecode = bindecoded';
       bindecode = bindecode(:);

     cont = cont +  (length(Binaire) - sum((Binaire == bindecode)));
      contt = contt +  (length(Binaire) - sum((Binaire == bindec)));
    end
    
 SyErr(s) = cont/(length(Binaire)*Npaquet);
  SyErrMP(s) = contt/(length(Binaire)*Npaquet);
end

save data2.mat




