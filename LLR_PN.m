%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LL=LLR_PN(y, hh, N0, Pvar, hMap, TYPE, SUM)%,Taw )

%%% calculates the LLRs in the phase-noise corrupted channel
% y - channel output (vector)
% hh- channel gain (vector the same size as y, or a scalar)
% N0 - variance of the complex Gaussian noise
% variance_PN - variance of the Gaussian phase noise, may be a vector
% hMap.X - constellation points
% hMap.pntk0 - indices to the symbols labeled with bit=0; a matrix, first column
%       corresponds to a first bit position;
%        number of columns may be smaller than log_2(M), LLRs are then calculated
%       only for the bits indicated by kk0 and kk1
% hMap.pntk1 - indices to the symbols labeled with bit=1
% TYPE - calculation method,
%       'LT' - linear transformation;
%       'BLT' - bilinear transformation;
%       'SP' -  approximate saddlepoint solution
%       'SP0' - Laplace transformation (saddlepoint method) with p=0
%       'SP1' - Laplace transformation (saddlepoint method) with p=delta
% SUM - combining the symbol metrics for LLR
%       'maxlog'  - consider two closest symbols
%       'sumexp'  - all symbols included
%

if nargin<7, SUM='sumexp'; end

variance_PN=Pvar;
% XX=hMap.Xsort(:);
XX=hMap.X(:);
M=length(XX);

kk0=hMap.pntk0;
kk1=hMap.pntk1;

[K,kM]=size(hMap.pntk0);
if K~=M/2, error('number of rows in kk0 must be equal to M/2'); end
LL=zeros(length(y),kM);

if length(variance_PN)==1
    variance_PN = variance_PN * ones( length(y),1 );
end
% Matrix-versions for simple operations
xref    = repmat( XX.',length(y),1 );
yy      = repmat( y,1,M );
% Taw      = repmat( Taw,1,M );
vVar_PN = repmat( variance_PN, 1, M );
switch TYPE
    case 'LT'
        y2 = yy; xref2 = xref;
        %%%% L_BLT_log ~=  log pdf(y|xref)
        LLR_log=-abs(y2-xref2).^2/N0;  %% Euclidean distance
        %%% add correction
        v2=N0+2*vVar_PN.*abs(xref2).^2;
        LLR_log=  LLR_log + imag(conj(xref2).*y2).^2*2.*(vVar_PN/N0)./v2 -0.5*log(v2);
%         LLR_log=LLR_log+ log(1+exp(log(Taw)-LLR_log));
        
    case 'BLT'
        y2 = (3*yy-xref)/2; xref2   = (yy+xref)/2;
        %%%% L_BLT_log ~=  log pdf(y|xref)
        LLR_log=-abs(y2-xref2).^2/N0;  %% Euclidean distance
        %%% add correction
        v2=N0+2*vVar_PN.*abs(xref2).^2;
        LLR_log=  LLR_log + imag(conj(xref2).*y2).^2*2.*(vVar_PN/N0)./v2 -0.5*log(v2);
        
    case 'SP0'
        y2 = yy; xref2 = xref;
        %%%% L_BLT_log ~=  log pdf(y|xref)
        LLR_log=-abs(y2-xref2).^2/N0;  %% Euclidean distance
        %%% add correction
        v2 = abs( N0+2*vVar_PN.*real(xref2.*conj(y2)) );
        LLR_log=  LLR_log + imag(conj(xref2).*y2).^2*2.*(vVar_PN/N0)./v2 -0.5*log(v2);
        
    case 'SP1'
        y2 = yy; xref2 = xref;
        %%%% L_BLT_log ~=  log pdf(y|xref)
        LLR_log=-abs(abs(y2)-abs(xref2)).^2/N0;  %% Euclidean distance
        
        delta = angle(y2./xref);
        v2 = abs( N0+2*vVar_PN.*abs(xref2.*y2) );
%         LLR_log=  LLR_log - 0.5*(delta).^2.*(2*abs(xref2.*y2))./v2 -0.5*log(v2) ;
        LLR_log=  LLR_log - 0.5.*delta.^2./vVar_PN + 0.5*delta.^2.*(N0./vVar_PN)./v2 -0.5*log(v2) ;
        
    case 'SP'
        y2 = yy; xref2 = xref;
        delta = angle(y2./xref); a=N0./(2*abs(xref.*y2).*vVar_PN);
        p_hat = delta./(1+a);
        
        LLR_log=-abs(y2-xref2.*exp(1j*p_hat)).^2/N0;  %% Euclidean distance
        
        v2 = abs( N0+2*vVar_PN.*abs(xref2.*y2).*cos(p_hat-delta) );
        LLR_log=  LLR_log - 0.5*log( v2) - 0.5.*p_hat.^2./vVar_PN;
        
    case 'SPit'
        y2 = yy; xref2 = xref;
        delta = angle(y2./xref); a=N0./(2*abs(xref.*y2).*vVar_PN);
        p_hat = delta./(1+a);
        for kk=1:2 %% Newton method to solve saddlepoint equation
            p_hat = p_hat - (sin(p_hat-delta) + p_hat.*a)./(cos(p_hat-delta) + a);
        end
        
        LLR_log=-abs(y2-xref2.*exp(1j*p_hat)).^2/N0;  %% Euclidean distance       
        v2 = abs( N0+2*vVar_PN.*abs(xref2.*y2).*cos(p_hat-delta) );
        LLR_log=  LLR_log - 0.5*log( v2) - 0.5.*p_hat.^2./vVar_PN;
    case 'EX'
        %% calculation using Gauss-Hermite numerical integration
        n=15;
        [p_GH, w_GH] = GaussHermite(n);
        
        gg_fun= @(p) exp(-abs( yy-xref.*exp( 1i*p*sqrt(vVar_PN*2) ) ).^2/N0)/(pi*sqrt(pi)*N0);
        
        LLR = zeros(size(yy));
        for nn=1:n
            LLR = LLR + gg_fun(p_GH(nn))*w_GH(nn);
        end
        LLR_log = log( LLR );
    otherwise
        error('unknown integral simplification');
end


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


