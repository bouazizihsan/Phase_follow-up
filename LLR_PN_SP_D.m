%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LL=LLR_PN_SP_D(y, N0,BCJR, hMap, SUM )

% q=[payload*S(2)-1 payload*S(2)].';
% q=q(:);

% variance =BCJR.Pvar(q,:);
variance =BCJR.Pvar;
a=variance(1:2:end,1);
b=variance(2:2:end,1);
d=variance(1:2:end,2);
c=variance(2:2:end,2);

if nargin<5, SUM='sumexp'; end

XX1=hMap(1).Xsort(:);
XX2=hMap(2).Xsort(:);
M=length(XX2);

kk0=hMap(2).pntk0;
kk1=hMap(2).pntk1;

[K,kM]=size(hMap(2).pntk0);
if K~=M/2, error('number of rows in kk0 must be equal to M/2'); end
LL=zeros(length(y),kM);

if length(d)==1
    d = d * ones( length(y),1 );
end
% Matrix-versions for simple operations
xref2    = repmat( XX2.',length(y),1 );
y2      = repmat( y(:,2),1,M );

xref1    = repmat( XX1.',length(y),1 );
y1      = repmat( y(:,1),1,M );
d = repmat( d, 1, M );
a      = repmat( a,1,M );
b      = repmat( b,1,M );
c      = repmat( c,1,M );

A1=(2*abs(xref1.*y1))/N0;
A2=(2*abs(xref2.*y2))/N0;

delta2 = angle(y2./xref2);delta1 = angle(y1./xref1);
[ind1,ind2]=find((c+b)==0);
p_hat2 = (2.*A2.*A1./((c+b).*(A1+a)).*delta2-delta1).*1./(2.*A1.*(A2+d)./((c+b).*(A1+a))- (c+b)./(2*A1) );
p_hat2(ind1,ind2)=1./(1+d(ind1,ind2)./A2(ind1,ind2)).*delta2(ind1,ind2);
LLR_log=-abs(y2-xref2.*exp(1j*p_hat2)).^2/N0;  %% Euclidean distance

v2 = abs( N0+2*d.*abs(xref2.*y2).*cos(p_hat2-delta2) );
LLR_log=  LLR_log - 0.5*log( v2) - 0.5.*p_hat2.^2./d;



for k=1:kM
    switch SUM
        case 'sumexp'
            LL(:,k)=jac_log( LLR_log(:,kk1(:,k)) ) - jac_log( LLR_log(:,kk0(:,k)) );
            
        case 'maxlog'
            LL(:,k)= max( LLR_log(:,kk1(:,k)), [],2 ) - max( LLR_log(:,kk0(:,k)), [],2);
            
        otherwise
            error('unknown metric combining');
    end
end
end


