function [y]= DECturbo(state,msgest,K,kk,n,L,B)

for i=1:L
for j=1:B

end
end


function [metric0,metric1]=mesure(state,msgest,L,n,K,B)
    metric=zeros(B,L);
    metric0=0;
    metric1=0;
    initialS=00;
for ii=1:L
     
            for m=1:n
   metric0=metric0+(msgest(ii,m) ~= 0)
   metric1=metric1+(msgest(ii,m) ~= 0)
            end
   metric(0,ii)=metric0;
   metric(1,ii)=metric1;
   
   if(ii>= (K-1))
            for m=1:n
   metric0=metric0+(msgest(ii,m) ~= 0)
   metric1=metric1+(msgest(ii,m) ~= 0)
            end
   metric(2,ii)=metric0;
   metric(3,ii)=metric0;
   
   end 
     
   
    end
end 

end