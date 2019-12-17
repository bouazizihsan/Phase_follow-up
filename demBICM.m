% Hsan ; = LLR_MQAM 
function  llr = demBICM(y,snr,M)
L = length(y);
mm=log2(M);
llr=zeros(L,mm);
for i=1:L
    lk=0;
    for k=1:mm % la position de kth bitb dans le symbole
        metric2=0;
        metric=0;
        
        for j = 1:M
            
            xx = de2bi(j-1,mm,'left-msb');
            if(xx(k))
                metric  = metric + exp(-snr*norm( y(i) - mapping(j-1,M)).^2);
            else
                metric2 = metric2 + exp(-snr*norm( y(i) - mapping(j-1,M)).^2);
            end
            
            
        end
        lk=lk - log(metric)+log(metric2);
        llr(i,k)=lk;
    end
   
end

end
% function  v = setting(M,k)
% x = qammod([0:M-1],M);
%
% end

function x = mapping(be,M)
x=0;
x = qammod(be,M);
end
