function LLR=LLR_Tikhonov1D(M,y_in,hMap1,hMap2,N0,BCJR1,SUM)
% LLR_Tikhonov calcule llR des bits pour une approximation d'un seul class

m=log2(M);
Ns=size(y_in,1);
LLR=zeros(Ns,m);
LLR_log=zeros(Ns,M);
Xref=[hMap1.Xsort hMap2.Xsort];
% Xref=Xref.';
% Xref=repmat(Xref,Ns,1);
% yy= repmat(y_in,1,M);
yy=y_in;
% Z1=(BCJR1.Z_alpha);
% Z2=(BCJR1.Z_beta);
Z=(BCJR1.Z_alpha_beta);
outbit0=hMap1.pntk0;
outbit1=hMap1.pntk1;

for i=1:M
% LLR_log(:,i)=sum(( abs(Z1 + Z2 + (yy.*repmat(conj(Xref(i,:)),Ns,1)*2./N0))-abs(Z1)- abs(Z2)-0.5*log( abs(Z1 + Z2 + (yy.*repmat(conj(Xref(i,:)),Ns,1)*2./N0 ))./(abs(Z1).*abs(Z2)) ) ),2);
LLR_log(:,i)=sum(( abs(Z + (yy.*repmat(conj(Xref(i,:)),Ns,1)*2./N0)) - abs(Z) - 0.5*log( abs(Z + (yy.*repmat(conj(Xref(i,:)),Ns,1)*2./N0 ))./(abs(Z)) ) ),2);
LLR_log(:,i) = sum(-abs(Xref(i,:)).^2/N0,2) + LLR_log(:,i);
end

% % LLR_log=[];
% % for j=1:size(Z,3)
% % for i=1:M
% % % LLR_log(:,i)=sum(( abs(Z1 + Z2 + (yy.*repmat(conj(Xref(i,:)),Ns,1)*2./N0))-abs(Z1)- abs(Z2)-0.5*log( abs(Z1 + Z2 + (yy.*repmat(conj(Xref(i,:)),Ns,1)*2./N0 ))./(abs(Z1).*abs(Z2)) ) ),2);
% % LLR_log=[LLR_log (sum(( abs(Z(:,:,j) + (yy.*repmat(conj(Xref(i,:)),Ns,1)*2./N0)) - abs(Z(:,:,j)) - 0.5*log( abs(Z(:,:,j) + (yy.*repmat(conj(Xref(i,:)),Ns,1)*2./N0 ))./(abs(Z(:,:,j))) ) ) - abs(Xref(i,:)).^2/N0,2))];
% % end
% % end

% aux= exp(-abs(x).^2 ./N0) .* (var2);

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
