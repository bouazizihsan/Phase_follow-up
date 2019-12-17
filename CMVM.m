function [Ztot]=CMVM(W,z,D) % [Ztot,corr]=CMVM(alph,z,D) %
% epsilon=0.01;
% alpha=alph./sum(alph);
% [ind val] = find(alpha>epsilon);
% z=z(ind,:);
% alpha=alpha(ind,:);
% cor=sum(alpha);
W=W./sum(W);
u=angle(sum(repmat(W,1,D).*(1-1./(2*abs(z) ) ).*exp(1j*angle(z)),1) );
k=1-(sum(repmat(W,1,D).*(1-1./(2*abs(z))).*cos(u-angle(z)),1)); 
% k=(sum(repmat(W,1,D)./abs(z),1));
k=1./(2*k);
% k=1./(k);
% if(1/(k) > Threshold)
%     K=0;
% end
Ztot=k.*exp(1j.*u);

end

