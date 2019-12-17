 clc;
clear all;
warning off;
snrdb = [-5:1:20];
snrLin = 10.^(snrdb/20);
N = 100;
x=randint(1,500);
%% 

%Mapping
for i=1:length(x)
    if(x(i)==1)
        xmod(i)=1;
    else
        xmod(i)=-1;
    end

end
dem = zeros(1,length(x));

for s= 1:length(snrdb)
    
    snr = snrLin(s);
    cont = 0;
    for h=1:N
      z  = randn(size(xmod))./ sqrt(snr);
      y  = xmod + z; 
      dem = (y>0); 
    
%     for i=1:length(y)
%         if(y(i)>=0)
%             dem(i)=1;
%         else
%             dem(i)=0;
%         end
%     end
   
     cont = cont +  (length(x)- sum((x==dem)));
    end
    
    bitr(s) = cont/(length(x)*N);% on divise par length(x) pour rendre bitr entre 0 et 1
      
end
%% 

semilogy(snrdb,bitr) 
xlabel('snr')
ylabel('berr')
grid on
save rel4.mat




