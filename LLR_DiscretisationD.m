function LLR=LLR_DiscretisationD(apo,hMap1,hMap2,SUM )

XX=[hMap1.Xsort hMap2.Xsort];
M=size(XX,1);
S=size(apo);
LLR_log=zeros(S(3),S(4));
k10=hMap1.pntk0;
k11=hMap1.pntk1;
% % k20=hMap2.pntk0;
% % k21=hMap2.pntk1;

[K1,kM]=size(hMap1.pntk0);
[K2,kM]=size(hMap2.pntk0);

LLR=zeros(S(3),kM);

for j=1:S(4)
LLR_log(:,j)=log(apo(1,1,:,j));
end
for k=1:kM
    switch SUM
        case 'sumexp'
            LLR(:,k)=jac_log( LLR_log(:,k11(:,k)) ) - jac_log( LLR_log(:,k10(:,k)) );          
        case 'maxlog'
            LLR(:,k)= max( LLR_log(:,k11(:,k)), [],2 ) - max( LLR_log(:,k10(:,k)), [],2);
        otherwise
            error('unknown metric combining');
    end
end
end
