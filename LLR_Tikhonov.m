function LLR=LLR_Tikhonov(hMap,M,y_in,x_in,N0,Z_alpha,Z_beta,fact_alpha,fact_beta,SUM)
% LLR_Tikhonov calcule llR des bits pour une approximation d'un class ou
% plus.
outbit0=hMap.pntk0;
outbit1=hMap.pntk1;

m=log2(M);
Ns=length(y_in);
LLR=zeros(Ns,m);
LLR_log=zeros(Ns,M);
Lmax=size(Z_alpha,2);
x=repmat(x_in,Ns,1);
Xref=repmat(x_in,Lmax,1);
Xref=Xref(:);
Xref=Xref.';
Xref=repmat(Xref,Ns,1);
yy= repmat(y_in,1,Lmax*M);
phi1=repmat(fact_alpha,1,M);
phi2=fact_beta;

Z1=repmat(Z_alpha,1,M);
Z2=Z_beta;

% var=zeros(Ns,Lmax*M);
%%somme sur les classes de alpha
% % for i=1:Lmax  
% %     var=var + phi1.*repmat(phi2(:,i),1,Lmax*M).*besseli(0, abs(Z1 + repmat(Z2(:,i),1,Lmax*M) + yy.*conj(Xref)*2./N0 ) )./ ( besseli(0, abs(Z1)).* besseli(0, abs(repmat(Z2(:,i),1,Lmax*M)) ));
% % end

var=[];
var2=[];
for j=1:Lmax:Lmax*(M-1)+1
  variable=[];
  for i=1:Lmax      
    var=[var log(phi1.*repmat(phi2(:,i),1,Lmax*M)) + ( abs(Z1 + repmat(Z2(:,i),1,Lmax*M) + yy.*conj(Xref)*2./N0 )-abs(Z1)- abs(repmat(Z2(:,i),1,Lmax*M))-0.5*log( abs(Z1 + repmat(Z2(:,i),1,Lmax*M) + yy.*conj(Xref)*2./N0 )./(abs(Z1).*abs(repmat(Z2(:,i),1,Lmax*M)))) )];    
  end
   if(Lmax>1)
      variable=[variable [var(:,j:j+(Lmax-1)) var(:,j+Lmax*M:j+Lmax*M+(Lmax-1))]];
      var2=[var2 jac_logH(variable)];
   else
      var2=[var2 [var(:,j:j+(Lmax-1))]];
      
   end
      
  
end



%%somme sur les classes de beta


LLR_log = (-abs(x).^2/N0) + (var2);

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
