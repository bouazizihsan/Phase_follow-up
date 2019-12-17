function LLR=LLR_Tikhonov1(M,y_in,hMap,N0,BCJR,SUM,payload)
% LLR_Tikhonov calcule llR des bits pour une approximation d'un seul class

m=log2(M);
Ns=length(y_in);
LLR=zeros(Ns,m);
LLR_log=zeros(Ns,M);
Xref=hMap.X(:).';
Xref=Xref(:);
Xref=Xref.';
Xref=repmat(Xref,Ns,1);
yy= repmat(y_in,1,M);

Z1=repmat(BCJR.Z_alpha(payload,:),1,M);
Z2=repmat(BCJR.Z_beta(payload,:),1,M);

outbit0=hMap.pntk0;
outbit1=hMap.pntk1;
% [outbit1,outbit0]=test(M);

% var=zeros(Ns,M);
%%somme sur les classes de alpha
% % for i=1:Lmax
% %     var=var + phi1.*repmat(phi2(:,i),1,Lmax*M).*besseli(0, abs(Z1 + repmat(Z2(:,i),1,Lmax*M) + yy.*conj(Xref)*2./N0 ) )./ ( besseli(0, abs(Z1)).* besseli(0, abs(repmat(Z2(:,i),1,Lmax*M)) ));
% % end

% k=abs(Z1) + abs(Z2) + abs(yy.*conj(Xref)*2./N0);
% ang=1./k.*(abs(Z1).*angle(Z1)+abs(Z2).*angle(Z2)+abs(yy.*conj(Xref)*2./N0).*angle(yy.*conj(Xref)*2./N0));
% Z=k.*exp(1j.*ang);

Z=(Z1) + (Z2) + (yy.*conj(Xref)*2./N0);
var=log(besseli(0,abs(Z),1))+abs(Z)-log(besseli(0,abs(Z1),1))-abs(Z1)-log(besseli(0,abs(Z2),1))-abs(Z2);
% var=( abs(Z)-abs(Z1)-abs(Z2)-0.5*log( abs(Z)./(abs(Z1).*abs(Z2)) ) );

LLR_log = (-abs(Xref).^2/N0) + var;

for k=1:m
    switch SUM
        case 'sumexp'
            LLR(:,k)=jac_logH( LLR_log(:,outbit1(:,k)) ) - jac_logH( LLR_log(:,outbit0(:,k)) );
            %LLR_exact(:,k)= sum(aux(:,outbit1(:,k)),2)- sum(aux(:,outbit0(:,k)),2);
        case 'maxlog'
            LLR(:,k) =( max(LLR_log(:,outbit1(:,k)).')-max(LLR_log(:,outbit0(:,k)).')).';
        otherwise
            error('unknown metric combining');
    end
end
end

% function [outbit1,outbit0]=test(n)
% m=log2(n);
% outbit1=zeros(n/2,m);
% outbit0=zeros(n/2,m);
% for k=1:m
%     b=bitget(0:n-1,m-k+1);% bitget(nombre decimal, bit), nombre decimal se transforme automatiquement en binaire, right-msb donc (m-k+1)
%     [row,col] = find(b);
%     outbit1(:,k)=col;
%     [row,col] = find(b==0);
%     outbit0(:,k)=col;
% end
% end
