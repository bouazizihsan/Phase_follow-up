function LL=LLR_SP1_P_Ref_XY(y, hh, N0, variance_PN, hMap, SUM, P_Ref )


if nargin<7, SUM='sumexp'; end

XX=hMap.Xsort(:);
M=length(XX);

kk0=hMap.pntk0;
kk1=hMap.pntk1;

[K,kM]=size(hMap.pntk0);
if K~=M/2, error('number of rows in kk0 must be equal to M/2'); end
LL=zeros(length(y),kM);
LL1=zeros(length(y),kM);

if length(variance_PN)==1
    variance_PN = variance_PN * ones( length(y),1 );
end
% Matrix-versions for simple operations
xref    = repmat( XX.',length(y),1 );
yy      = repmat( y,1,M );
vVar_PN = repmat( variance_PN, 1, 1, M );
P_Reff=repmat(P_Ref,1,M);
% P_Reff=angle(yy./xref);
%case 'SP1'
y2 = conj(yy).*exp(1j*P_Reff); xref2 = xref;
LLR_log=-0.5*log(pi*N0*diag(abs(xref2).^2))-(abs(real(y2.*xref2)-(abs(xref2).^2))).'*inv(N0*diag(abs(xref2).^2))*(abs(real(y2.*xref2)-(abs(xref2).^2)));  %% Euclidean distance

delta =P_Reff - imag(y2./conj(xref));
v2 =  (N0./2)*inv(diag(abs(xref2).^2))+vVar_PN ;
LLR_log=LLR_log-0.5*log(v2)-(delta).^2./(2.*v2) ;

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


