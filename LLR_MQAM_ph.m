function LLR=LLR_MQAM_ph(M,y_in,x_in,N0,sig,LLR_METHOD,SUM)
m=log2(M);
s=size(y_in);
LLR=zeros(s(1).*s(2),m);

if length(sig)==1
    sig = sig * ones( length(y_in),1 );
end
LLR_log=zeros(length(y_in),M);
[outbit1,outbit0]=test(M);
xref= repmat(x_in,s(1).*s(2),1);
yy= repmat(y_in,1,M);
sigw2 = repmat( sig, 1, M );
p0= angle(yy./xref);
a=N0./(2.*abs(yy).*abs(xref).*sigw2);
p_tild=p0./(1+a);

switch LLR_METHOD
    case 'LT'
        LLR_log =  ((-abs(yy -xref).^2 ./N0)  +  (imag(conj(xref).*yy ).^2.* 2.*sigw2./N0)./ (2.*sigw2.*abs(xref).^2 + N0))  -(0.5).*log( N0 + 2.*sigw2.*abs(xref).^2 ) ;
    case 'BLT'
        LLR_log= ((-abs(yy -xref).^2 .*(1./N0))   +  (imag(conj(xref).*yy ).^2.* 2.*sigw2.*(1./N0)./ (0.5.*sigw2.*abs(xref+yy ).^2+1./(1./N0)))  -(0.5).*log((1./((1./N0)))+ 0.5.*sigw2.*abs(xref+ yy ).^2 )) ;
    case 'SP'
        LLR_log=   ((-abs(yy -xref.*exp(1j*p_tild)).^2 .*(1./N0)) + (-p_tild.^2./(2.*sigw2)) -(0.5).*log( abs((2.*abs(yy ).*abs(xref).* cos(p_tild-p0).*sigw2 + 1./((1./N0)))./(sigw2./((1./N0)))  ) ));
        
    case 'SP0'
        LLR_log=  ((-abs(yy -xref).^2./N0)   + (imag(conj(xref).*yy ).^2.* 2.*(sigw2/N0)./ (2.*sigw2.*real(conj(yy).* xref )+ N0))   -(0.5).*log((1./((1./N0)))+ 2.*sigw2.*abs(real(conj(xref).* yy ) )));
    
    case 'SP1'
        LLR_log=   ((-abs(abs(yy )-abs(xref)).^2 .*(1./N0)) + (-p0.^2./(2.*sigw2)) -(0.5).*log((2.*abs(yy ).*abs(xref).*sigw2 + 1./((1./N0)))./(sigw2./((1./N0)))  ) + (p0.^2./(sigw2 .*(1./N0))./(2.*(2.*abs(yy ).*abs(xref).*sigw2+ 1./((1./N0))) ) ));
    
    otherwise
        error('unknown integral simplification');
end
for k=1:m
    switch SUM
        case 'sumexp'
            LLR(:,k)=jac_logH( LLR_log(:,outbit1(:,k)) ) - jac_logH( LLR_log(:,outbit0(:,k)) );
        case 'maxlog'
            LLR(:,k) =( max(LLR_log(:,outbit1(:,k)).')-max(LLR_log(:,outbit0(:,k)).')).';
        otherwise
            error('unknown metric combining');
    end
end
end

function [outbit1,outbit0]=test(n)
m=log2(n);
outbit1=zeros(n/2,m);
outbit0=zeros(n/2,m);
for k=1:m
    b=bitget(0:n-1,m-k+1);% bitget(nombre decimal, bit), nombre decimal se transforme automatiquement en binaire, right-msb donc (m-k+1)
    [row,col] = find(b);
    outbit1(:,k)=col;
    [row,col] = find(b==0);
    outbit0(:,k)=col;
end
end

