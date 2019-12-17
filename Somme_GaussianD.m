function [moy, var]=Somme_GaussianD(m,sigma2,phi)
phi=phi./sum(phi);
s=size(m);
D=s(2);
sumphi=sum(phi);
moy=sum(repmat(phi,1,D).*m,1)./sum(phi,1);
% for j=1:s(1)
%     sig(j,:)=diag(sigma2(1+(j-1)*D:j*D,:))';
% end
% moy=angle( sum( repmat(phi,1,D).*( besseli(1,1./sig,1)./besseli(0,1./sig,1).*exp(1j.*m) ),1 ) );%angle( ( phi'*( (1-sigma/2).*exp(1j.*m) ) )./sumWeight );

variable=zeros(D,D);
for j=1:s(1)
    variable=variable + phi(j).* (  ((m(j,:)-moy)).'*((m(j,:)-moy)) + sigma2(1+(j-1)*D:j*D,:) );
%     variable=variable + phi(j).* (  (Operator(m(j,:)-moy)).'*(Operator(m(j,:)-moy)) + sigma2(1+(j-1)*D:j*D,:) );
end
var=variable./(sumphi);
end