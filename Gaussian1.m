function [U_min, sigma2_min]=Gaussian1(m,sigma,phi)
Phase_limit=pi;
L=length(m);
phi=phi./(sum(phi));
[m,fs]= sort(m);
sigma=sigma(fs);
phi=phi(fs);

m2=[m;m+2.*Phase_limit];
sigma2=[sigma;sigma];

phi2=[phi;phi];

moy=zeros(size(m));
var=zeros(size(m));



for i=1:L
    
    moy(i)=sum(phi2(i:L+i-1).*m2(i:L+i-1))./(sum(phi2(i:L+i-1)));
    
    var(i)=sum(phi2(i:L+i-1).*(m2(i:L+i-1).^2+ sigma2(i:L+i-1)))./(sum(phi2(i:L+i-1))) -  moy(i).^2;
    
    %     if((i==1))
    %       figure(1)
    %       plot(pp,pdf_Gauss(pp,moy(1),var(1))./sum(pdf_Gauss(pp,moy(1),var(1))) , 'k');
    %       hold on
    %    end
    
end
[sigma2_min,indx]=min(var);
U_min=moy(indx);

U_min(U_min>Phase_limit)=U_min(U_min>Phase_limit)-2*Phase_limit;

U_min(U_min<-Phase_limit)=U_min(U_min<-Phase_limit)+2*Phase_limit;



% if(abs(U_min-Phase_noise) < abs(moy(1)-Phase_noise))
%   stp=1
% end
% [abs(U_min-Phase_noise) abs(moy(1)-Phase_noise)]
%
% figure(1)
% hold on
% aaa=phi'*pdf_Gauss(repmat(pp,1,length(m)), repmat(m',length(pp),1 ), repmat(sigma',length(pp),1 ) ).';
% plot(pp,aaa./sum(aaa) , 'b');
%
% plot(pp,pdf_Gauss(pp,U_min,sigma2_min)./sum(pdf_Gauss(pp,U_min,sigma2_min)) , 'g--');
% hold on
% plot(Phase_noise,0,'r*')
% hold off


end