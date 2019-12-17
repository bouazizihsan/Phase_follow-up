% pour le codage.m
function l=LLR_MQAM(IM,y_in,x_in,SNR)

K=log2(IM);

s=size(y_in);

l=zeros(s(1)*s(2),K);

for i =1:s(1)*s(2)
    
    for k=1:K
        
        [outbit1,outbit0]=test(x_in,k,IM); % test sur le k ieme bit de x
        
        num=0;
        
        denum=0;
        
        for j=1:IM
            
            if abs(outbit1(j))
                
                num = num + exp(-SNR*abs(y_in(i)-outbit1(j)).^2); % calcul des metrics
                
            else
                denum= denum + exp(-SNR*abs(y_in(i)-outbit0(j)).^2);
                
            end
            
        end
        
        l(i,k) = log(num/denum);
        
    end
    
end

end
function [outbit1,outbit0]=test(x,k,n)
outbit1=zeros(1,n);
outbit0=zeros(1,n);
K=log2(n);
for m=1:n
    b=bitget(m-1,K-k+1);% bitget(nombre decimal, bit), nombre decimal se transforme automatiquement en binaire
    if b
        outbit1(m)=x(m);
    else
        outbit0(m)=x(m);
    end
end
end

