% pour le codage.m
function [LT,BLT,SP,SP0,SP1]=LLR_MQAMphase(IM,y_in,x_in,snr,sigma_p2)

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
        
        num=0;
        numBLT=0;
        denum=0;
        denumBLT=0;
        numSP=0;
        denumSP=0;
        numSP0=0;
        denumSP0=0;
        numSP1=0;
        denumSP1=0;
        
        for jj=1:IM
            
            if abs(outbit1(jj))
                p0= angle(y_in(i)./outbit1(jj));
                
                a=1./(snr *2*abs(y_in(i)).*abs(outbit1(jj)).*sigma_p2);
                p_tild =p0./(1+a);
                %num =[num  -2*snr*abs(y_in(i)-outbit1(j)).^2  +  ((2*sigma_p*2*snr)./ (2*sigma_p*abs(outbit1(j)).^2  +  1./(2*snr))) *imag(conj(outbit1(j))*y_in(i)).^2 - 0.5*log((1./(2*snr))+(2*sigma_p*abs(outbit1(j)).^2))]; % calcul des metrics
                num= num + exp(-abs(y_in(i)-outbit1(jj)).^2 *snr)   *  exp(imag(conj(outbit1(jj))*y_in(i)).^2* 2*sigma_p2*snr./ (2*sigma_p2*abs(outbit1(jj)).^2+1./snr))  ./ ((1./(snr))+ 2*sigma_p2*abs(outbit1(jj)).^2 ).^(0.5);
                
                
                %numBLT=[numBLT -2*snr*abs(y_in(i)-outbit1(j)).^2  +  ((2*sigma_p*2*snr)./ (2*sigma_p*real( conj(outbit1(j)) .* y_in(i))   +  1./(2*snr))) *imag(conj(outbit1(j))*y_in(i)).^2 - 0.5*log((1./(2*snr))+(2*sigma_p*real(conj(outbit0(j))* y_in(i) ))]; % calcul des metrics
                numBLT= numBLT + exp(-abs(y_in(i)-outbit1(jj)).^2 *snr)   *  exp(imag(conj(outbit1(jj))*y_in(i)).^2* 2*sigma_p2*snr./ (0.5*sigma_p2*abs(outbit1(jj)+ y_in(i)).^2+1./snr))  ./ ((1./(snr))+ 0.5*sigma_p2*abs(outbit1(jj)+ y_in(i)).^2 ).^(0.5);
                
                numSP= numSP + exp(-abs(y_in(i)-outbit1(jj)* exp(1j*p_tild)).^2 *snr) *exp(-p_tild.^2./(2*sigma_p2)) ./ sqrt( abs((2*abs(y_in(i))*abs(outbit1(jj))* cos(p_tild-p0)*sigma_p2 + 1./(snr)./(sigma_p2./(snr)))  ) );
                
                
                numSP0= numSP0 + exp(-abs(y_in(i)-outbit1(jj)).^2 *snr)   *  exp(imag(conj(outbit1(jj))*y_in(i)).^2* 2*sigma_p2*snr./ (2*sigma_p2*real(conj(outbit1(jj)).* y_in(i))+1./snr))  ./ ((1./(snr))+ 2*sigma_p2*real(conj(outbit1(jj)).* y_in(i)) ).^(0.5);
                
                
                numSP1= numSP1 + exp(-abs(abs(y_in(i))-abs(outbit1(jj))).^2 *snr) *exp(-p0.^2./(2*sigma_p2))  ./sqrt((2*abs(y_in(i))*abs(outbit1(jj))*sigma_p2 + 1./(snr))./(sigma_p2./(snr))) *exp(p0.^2./(sigma_p2 *snr)./(2*(2*abs(y_in(i))*abs(outbit1(jj))*sigma_p2+ 1./(snr)) ) );
                
            else
                p0= angle(y_in(i)./outbit0(jj));
                a=1./(snr *2*abs(y_in(i)).*abs(outbit0(jj)).*sigma_p2);
                p_tild =p0./(1+a);
                %denum= [denum -2*snr*abs(y_in(i)-outbit0(j)).^2+((2*sigma_p*2*snr)./(2*sigma_p*abs(outbit0(j)).^2 + 1./(2*snr)))*imag(conj(outbit0(j))*y_in(i)).^2 -0.5*log( (1./(2*snr)) +(2*sigma_p*abs(outbit0(j)).^2))];
                denum= denum + exp(-abs(y_in(i)-outbit0(jj)).^2 *snr)   *  exp(imag(conj(outbit0(jj))*y_in(i)).^2* 2*sigma_p2*snr./ (2*sigma_p2*abs(outbit0(jj)).^2+1./snr))  ./ ((1./(snr))+ 2*sigma_p2*abs(outbit0(jj)).^2 ).^(0.5);
                
                
                %denumBLT=[denumBLT -2*snr*abs(y_in(i)-outbit0(j)).^2+((2*sigma_p*2*snr)./(2*sigma_p*real(conj(outbit0(j)).* y_in(i)) )   + 1./(2*snr)))*imag(conj(outbit0(j))*y_in(i)).^2 -0.5*log( (1./(2*snr)) +(2*sigma_p*real(conj(outbit0(j))* y_in(i) )))];
                
                denumBLT= denumBLT + exp(-abs(y_in(i)-outbit0(jj)).^2 *snr)   *  exp(imag(conj(outbit0(jj))*y_in(i)).^2* 2*sigma_p2*snr./ (0.5*sigma_p2*abs(outbit0(jj)+ y_in(i)).^2+1./snr)) ./ ((1./(snr))+ 0.5*sigma_p2*abs(outbit0(jj)+ y_in(i)).^2 ).^(0.5);
                
                denumSP=denumSP + exp(-abs(y_in(i)-outbit0(jj)*exp(1j*p_tild)).^2 *snr) *exp(-p_tild.^2./(2*sigma_p2)) ./ sqrt( abs((2*abs(y_in(i))*abs(outbit0(jj))* cos(p_tild-p0)*sigma_p2 + 1./(snr))./(sigma_p2./(snr))  ) );
                
                denumSP0= denumSP0 + exp(-abs(y_in(i)-outbit0(jj)).^2 *snr)   *  exp(imag(conj(outbit0(jj))*y_in(i)).^2* 2*sigma_p2*snr./ (2*sigma_p2*real(conj(outbit0(jj)).* y_in(i))+1./snr)) ./ ((1./(snr))+ 2*sigma_p2*real(conj(outbit0(jj)).* y_in(i)) ).^(0.5);
                
                denumSP1= denumSP1 + exp(-abs(abs(y_in(i))-abs(outbit0(jj))).^2 *snr) *exp(-p0.^2./(2*sigma_p2))  ./sqrt((2*abs(y_in(i))*abs(outbit0(jj))*sigma_p2 + 1./(snr))./(sigma_p2./(snr)) ) *exp(p0.^2./(sigma_p2 *snr)./(2*(2*abs(y_in(i))*abs(outbit0(jj))*sigma_p2+ 1./(snr)) ) );
                
            end
            
        end
       
        LT(i,k) = log(num./denum);
        BLT(i,k)= log(numBLT./denumBLT);
        SP(i,k)= log(numSP./denumSP);
        SP0(i,k)= log(numSP0./denumSP0);
        SP1(i,k)= log(numSP1./denumSP1);
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

