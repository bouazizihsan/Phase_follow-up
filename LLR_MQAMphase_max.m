% pour le codage.m
function [LT,BLT,SP,SP0,SP1]=LLR_MQAMphase_max(IM,y_in,x_in,snr,sigma_p)

K=log2(IM);

s=size(y_in);
LT=zeros(s(1)*s(2),K);
BLT=zeros(s(1)*s(2),K);
SP=zeros(s(1)*s(2),K);
SP0=zeros(s(1)*s(2),K);
SP1=zeros(s(1)*s(2),K);

for i =1:s(1)*s(2)
    
    for k=1:K
        
        [outbit1,outbit0]=test(x_in,k,IM); % test sur le k ieme bit de x
        
        num=[];
        numBLT=[];
        denum=[];
        denumBLT=[];
        numSP=[];
        denumSP=[];
        numSP0=[];
        denumSP0=[];
        numSP1=[];
        denumSP1=[];
        
        for j=1:IM
            
            if abs(outbit1(j))
                p0= angle(y_in(i)./outbit1(j));
                
                a=1./(snr *2*abs(y_in(i)).*abs(outbit1(j)).*sigma_p.^2);
                p_tild =p0./(1+a);
                num =[num,  ((-abs(y_in(i)-outbit1(j)).^2 *snr)  +  (imag(conj(outbit1(j))*y_in(i)).^2* 2*sigma_p.^2*snr./ (2*sigma_p.^2*abs(outbit1(j)).^2+1./snr))  -(0.5)*log( ((1./(snr))+ 2*sigma_p.^2*abs(outbit1(j)).^2 )))];
                
                numBLT=[numBLT, ((-abs(y_in(i)-outbit1(j)).^2 *snr)   +  (imag(conj(outbit1(j))*y_in(i)).^2* 2*sigma_p.^2*snr./ (0.5*sigma_p.^2*abs(outbit1(j)+y_in(i)).^2+1./snr))  -(0.5)*log((1./(snr))+ 0.5*sigma_p.^2*abs(outbit1(j)+ y_in(i)).^2 )) ];
                
                numSP= [numSP,  ((-abs(y_in(i)-outbit1(j)*exp(1j*p_tild)).^2 *snr) + (-p_tild.^2./(2*sigma_p.^2)) -(0.5)*log( abs(2*abs(y_in(i))*abs(outbit1(j))* cos(p_tild-p0)*sigma_p.^2 + 1./(snr)  ) ))];
                
                
                numSP0= [numSP0, ((-abs(y_in(i)-outbit1(j)).^2 *snr)   + (imag(conj(outbit1(j))*y_in(i)).^2* 2*sigma_p.^2*snr./ (2*sigma_p.^2*real(conj(outbit1(j)).* y_in(i))+1./snr))   -(0.5)*log((1./(snr))+ 2*sigma_p.^2*real(conj(outbit1(j)).* y_in(i)) ))];
                
                
                numSP1= [numSP1,  ((-abs(abs(y_in(i))-abs(outbit1(j))).^2 *snr) + (-p0.^2./(2*sigma_p.^2)) -(0.5)*log(2*abs(y_in(i))*abs(outbit1(j))*sigma_p.^2 + 1./(snr)  ) + (p0.^2*sigma_p.^2 *snr./(2*(2*abs(y_in(i))*abs(outbit1(j))*sigma_p.^2+ 1./(snr)) ) ))];
                
            else
                p0= angle(y_in(i)./outbit0(j));
                a=1./(snr *2*abs(y_in(i)).*abs(outbit0(j)).*sigma_p.^2);
                p_tild =p0./(1+a);
                denum =[denum,  ((-abs(y_in(i)-outbit0(j)).^2 *snr)  +  (imag(conj(outbit0(j))*y_in(i)).^2* 2*sigma_p.^2*snr./ (2*sigma_p.^2*abs(outbit0(j)).^2+1./snr))  -(0.5)*log( ((1./(snr))+ 2*sigma_p.^2*abs(outbit0(j)).^2 )))];
                
                denumBLT=[denumBLT, ((-abs(y_in(i)-outbit0(j)).^2 *snr)   +  (imag(conj(outbit0(j))*y_in(i)).^2* 2*sigma_p.^2*snr./ (0.5*sigma_p.^2*abs(outbit0(j)+ y_in(i)).^2+1./snr))  -(0.5)*log((1./(snr))+ 0.5*sigma_p.^2*abs(outbit0(j)+ y_in(i)).^2 ))];
                
                denumSP= [denumSP,  ((-abs(y_in(i)-outbit0(j)*exp(1j*p_tild)).^2 *snr) + (-p_tild.^2./(2*sigma_p.^2)) -(0.5)*log( abs(2*abs(y_in(i))*abs(outbit0(j))* cos(p_tild-p0)*sigma_p.^2 + 1./(snr)  ) ))];
                
                
                denumSP0= [denumSP0, ( (-abs(y_in(i)-outbit0(j)).^2 *snr)   + (imag(conj(outbit0(j))*y_in(i)).^2* 2*sigma_p.^2*snr./ (2*sigma_p.^2*real(conj(outbit0(j)).* y_in(i))+1./snr))   -(0.5)*log((1./(snr))+ 2*sigma_p.^2*real(conj(outbit0(j)).* y_in(i)) ))];
                
                
                denumSP1= [denumSP1, ( (-abs(abs(y_in(i))-abs(outbit0(j))).^2 *snr) + (-p0.^2./(2*sigma_p.^2)) -(0.5)*log(2*abs(y_in(i))*abs(outbit0(j))*sigma_p.^2 + 1./(snr)  ) + (p0.^2*sigma_p.^2 *snr./(2*(2*abs(y_in(i))*abs(outbit0(j))*sigma_p.^2+ 1./(snr)) ) ))];
                
            end
            
        end
        num1=max(num);
        denum0= max (denum);
        numBLT1=max(numBLT);
        denumBLT0= max (denumBLT);
        
        numSP_1=max(numSP);
        denumSP_0= max (denumSP);
        
        numSP01=max(numSP0);
        denumSP00= max (denumSP0);
        
        numSP11=max(numSP1);
        denumSP10= max (denumSP1);
        
        LT(i,k) = (num1-denum0);
        BLT(i,k)= (numBLT1-denumBLT0);
        
        SP(i,k)= numSP_1-denumSP_0;
        SP0(i,k)= (numSP01-denumSP00);
        SP1(i,k)= (numSP11-denumSP10);
    end
    
end

end

function [outbit1,outbit0]=test(x,k,n)
outbit1=zeros(1,n);
outbit0=zeros(1,n);
K=log2(n);
for m=1:n
    b=bitget(m-1,K-k+1);% bitget(nombre decimal, bit), nombre decimal se transforme automatiquement en binaire, right-msb donc (m-k+1)
    if b
        outbit1(m)=x(m);
    else
        outbit0(m)=x(m);
    end
end
end

