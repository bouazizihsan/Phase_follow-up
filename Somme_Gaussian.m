function [moy, var]=Somme_Gaussian(m,sigma2,phi)
% m,sigma,phi vectors (M,1)

sumWeight=sum(phi);
phi=phi./sumWeight;
% moy=(phi'*m)./sumWeight;
% var=(phi'*(m.^2 + sigma2))./sumWeight -  moy.^2;

moy=angle( ( phi'*( besseli(1,1./sigma2,1)./besseli(0,1./sigma2,1).*exp(1j.*m) ) ) );%angle( ( phi'*( (1-sigma/2).*exp(1j.*m) ) )./sumWeight );
var=(phi'*(sigma2+(Operator(m-moy)).^2))./sumWeight ;

% % % var=2*(1-phi'*( (1-sigma2./2).*cos(moy-m) ) ) ; 

end