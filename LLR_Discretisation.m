function LLR=LLR_Discretisation(apo,hMap1,SUM)
%% size(apo)= (Ns,M)
XX=hMap1.X;
M=size(XX,1);
S=size(apo);
LLR_log=zeros(S(1),S(2));
k10=hMap1.pntk0;
k11=hMap1.pntk1;
[K1,kM]=size(hMap1.pntk0); 

LLR=zeros(S(1),kM);

% LLR_log=log(apo);
% 
% for k=1:kM
%     switch SUM
%         case 'sumexp'
%             LLR(:,k)=jac_log( LLR_log(:,k11(:,k)) ) - jac_log( LLR_log(:,k10(:,k)) );          
%         case 'maxlog'
%             LLR(:,k)= max( LLR_log(:,k11(:,k)), [],2 ) - max( LLR_log(:,k10(:,k)), [],2);
%         otherwise
%             error('unknown metric combining');
%     end
% end
for k=1:kM   
    LLR(:,k)=(log(sum(apo(:,k11(:,k)),2)) - log(sum(apo(:,k10(:,k)),2)));
end
end
