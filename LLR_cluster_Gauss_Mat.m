% function LLR=LLR_cluster_Gauss_Mat(M,y_in,x_in,N0,parametrs1,parametrs2,SUM,TYPE)
function LLR=LLR_cluster_Gauss_Mat(M,y_in,x_in,N0,G1,G2,SUM,TYPE)
m=log2(M);
Ns=length(y_in);
LLR=zeros(Ns,m);
% Lmax=size(parametrs1,2)/3;
Lmax=size(G1.MG1,2);
LLR_log=zeros(Ns,M*Lmax*Lmax);
Coef= zeros(Ns,Lmax^2);

[outbit1,outbit0]=test(M);

moy1=G1.MG1;
moy2=G2.MG2;
var1=G1.VG1;
var2=G2.VG2;
phi1=G1.PoidG1;
phi2=G2.PoidG2;
% moy1=parametrs1(:,1:Lmax);
% moy2=parametrs2(:,1:Lmax);
% var1=parametrs1(:,Lmax+1:2*Lmax);
% var2=parametrs2(:,Lmax+1:2*Lmax);
% phi1=parametrs1(:,2*Lmax+1:3*Lmax);
% phi2=parametrs2(:,2*Lmax+1:3*Lmax);

y= repmat(y_in,1,Lmax^2);


c=phi1';
c=repmat(c,1,Lmax);
c=c.';
c=c(:);
c=c.';
phi1=reshape(c,Ns,length(c)/Ns);

c=var1';
c=repmat(c,1,Lmax);
c=c.';
c=c(:);
c=c.';
var1=reshape(c,Ns,length(c)/Ns);

c=moy1';
c=repmat(c,1,Lmax);
c=c.';
c=c(:);
c=c.';
moy1=reshape(c,Ns,length(c)/Ns);

phi2=repmat(phi2,1,Lmax);
moy2=repmat(moy2,1,Lmax);
var2=repmat(var2,1,Lmax);

v=var2+var1;
var_in= (var2.*var1)./v;
u_in= (moy2.*var1 + moy1.*var2)./v;

for i=1:Lmax^2
    a=find(var1(:,i)==inf);
    b=find(var2(:,i)==inf);
    var_in(a,i)=var2(a,i);   %%% INF case of pilots or mixture no classified 
    var_in(b,i)=var1(b,i);
    u_in(a,i)=moy2(a,i);
    u_in(b,i)=moy1(b,i);

    c=find((var1(:,i)==0) + (var2(:,i)==0));
    var_in(c,i)=0;
    u_in(c,i)=0;
end


Coef=phi1.*phi2;% size= (Ns * L^2)
y= y.*exp(-1j.*u_in);

yy= repmat(y,1,M);
xref=repmat(x_in,Lmax^2,1);
xref=xref(:);
xref=xref.';
xref=repmat(xref,Ns,1);
vVar_PN=repmat(var_in,1,M);
%Coef=Coef./repmat(sum(Coef),Ns,1);

switch TYPE
    case 'LT'
        y2 = yy; xref2 = xref;
        %%%% L_BLT_log ~=  log pdf(y|xref)
        LLR_log=-abs(y2-xref2).^2/N0;  %% Euclidean distance
        %%% add correction
        v2=N0+2*vVar_PN.*abs(xref2).^2;
        LLR_log=  LLR_log + imag(conj(xref2).*y2).^2*2.*(vVar_PN/N0)./v2 -0.5*log(v2);
        
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


for k=1:m
    vect1=[];
    vect0=[];
    for f=1:size(outbit1,1)
        vect1=[vect1 [(outbit1(f,k)-1)*Lmax^2+1:outbit1(f,k)*Lmax^2]];
        vect0=[vect0 [(outbit0(f,k)-1)*Lmax^2+1:outbit0(f,k)*Lmax^2]];
    end
    switch SUM
        case 'sumexp'
            LLR(:,k)=jac_log_Nclass2( LLR_log(:,vect1 ),Coef,Lmax^2,M) - jac_log_Nclass2( LLR_log(:,vect0 ),Coef,Lmax^2,M);
        case 'maxlog'
            LLR(:,k) = (max(LLR_log(:,vect1).')- max(LLR_log(:,vect0).')).';
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
