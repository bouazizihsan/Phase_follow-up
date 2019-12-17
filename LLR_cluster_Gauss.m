function [LLR]=LLR_cluster_Gauss(M,y_in,x_in,N0,L_alpha,L_beta,parametrs1,parametrs2,SUM,TYPE)

m=log2(M);
Ns=length(y_in);
LLR=zeros(Ns,m);
%Lmax=size(parametrs1,2)./3;
L=size(parametrs1,2)./3;
LLR_log=zeros(Ns,M*L.^2);
Coef= zeros(Ns,L.^2);
Coef_INF= zeros(Ns,(L-1).^2);
%Xref=repmat(x_in,s,1);
%yy= repmat(y_in,1,M);
[outbit1,outbit0]=test(M);
pdf_Gauss   = @(t,m,s2) exp( -abs(t-m ).^2./(2*s2) )./sqrt(2*pi*s2); % pdf(t|m), normpdf(x,mu,sigma.^2)

for n=1:Ns
    aux=find(parametrs1(n,L+1:L+L_alpha(n))==Inf);
    
    fact=[];
    fact_INF=[];
    var_in=[];
    y=[]; 
    test_INF=0;
    
     for i=1:L_beta(n) %somme sur beta
        if((parametrs2(n,L+i)==Inf) && norm(aux))
            test_INF=aux;
            v=parametrs2(n,L+i)+parametrs1(n,L+1:L+L_alpha(n));
            u_jn= (parametrs1(n,1:L_alpha(n)).*parametrs2(n,L+i) + parametrs1(n,L+1:L+L_alpha(n)).*parametrs2(n,i))./v;
            var=(parametrs2(n,L+i).*parametrs1(n,L+1:L+L_alpha(n)))./v;
            u_jn(aux)=[];
            var(aux)=[];
            var_in=[var_in var];
        elseif(parametrs2(n,L+i)==Inf)
            u_jn=parametrs1(n,1:L_alpha(n));
            var=parametrs1(n,L+1:L+L_alpha(n));
            var_in=[var_in var];
        elseif(norm(aux))
            v=parametrs2(n,L+i)+parametrs1(n,L+1:L+L_alpha(n));
            u_jn= (parametrs1(n,1:L_alpha(n)).*parametrs2(n,L+i) + parametrs1(n,L+1:L+L_alpha(n)).*parametrs2(n,i))./v;
            var=(parametrs2(n,L+i).*parametrs1(n,L+1:L+L_alpha(n)))./v;
            u_jn(aux)=parametrs2(n,i);
            var(aux)=parametrs2(n,L+i);
            var_in=[var_in var];
        else
            v=parametrs2(n,L+i)+parametrs1(n,L+1:L+L_alpha(n));
            u_jn= (parametrs1(n,1:L_alpha(n)).*parametrs2(n,L+i) + parametrs1(n,L+1:L+L_alpha(n)).*parametrs2(n,i))./v;
            var_in=[var_in ((parametrs2(n,L+i).*parametrs1(n,L+1:L+L_alpha(n)))./v)];
        end
        
        fact=[fact parametrs2(n,2*L+i).*parametrs1(n,2*L+1:2*L+L_alpha(n)).*pdf_Gauss(parametrs2(n,i),parametrs1(n,1:L_alpha(n)),parametrs2(n,L+i)+parametrs1(n,L+1:L+L_alpha(n)))];
        if(test_INF)  
           fact_INF=[fact_INF fact(test_INF)];
           fact(test_INF)=[];
        end
        
        y=[y y_in(n).*exp(-1j.*u_jn)]; 
        
     end
    y=y(:);%L_alpha*L_bata
    var_in=var_in(:);
    ll=length(y);
    xref=repmat(x_in,ll,1);
%     Xref=Xref(:);
%     Xref=Xref.';
    yy= repmat(y,1,M);
    vVar_PN=repmat(var_in,1,M);
    %fact=fact./sum(fact);
    Coef_INF(n,1:length(fact_INF)) = fact_INF;
    Coef(n,1:ll) = fact;
    switch TYPE
        case 'LT'
            y2 = yy; xref2 = xref;
            
            %%%% L_BLT_log ~=  log pdf(y|xref)
            LLR_lo=-abs(y2-xref2).^2/N0;  %% Euclidean distance
            LLR_lo=[LLR_lo;zeros(L-ll,M)]; %size(LLR_lo,2)==M
            LLR_lo= LLR_lo(:);
            LLR_lo=LLR_lo';
            %%% add correction
            v2=N0+2*vVar_PN.*abs(xref2).^2;
            VV=imag(conj(xref2).*y2).^2*2.*(vVar_PN/N0)./v2 -0.5*log(v2);
            VV=[VV;zeros(L-ll,M)];
            VV=VV(:);
            VV=VV';
            LLR_log(n,:)=  LLR_lo + VV;
            
        case 'BLT'
            y2 = (3*yy-xref)/2; xref2   = (yy+xref)/2;
            %%%% L_BLT_log ~=  log pdf(y|xref)
            LLR_lo=-abs(y2-xref2).^2/N0;  %% Euclidean distance
            LLR_lo=[LLR_lo;zeros(L.^2-ll,M)];
            LLR_lo= LLR_lo(:);
            LLR_lo=LLR_lo.';
            %%% add correction
            v2=N0+2*vVar_PN.*abs(xref2).^2;
            VV=imag(conj(xref2).*y2).^2*2.*(vVar_PN/N0)./v2 -0.5*log(v2);
            VV=[VV;zeros(L.^2-ll,M)];
            VV=VV(:);
            VV=VV.';
            LLR_log(n,:)=  LLR_lo + VV;
            
        case 'SP0'
            y2 = yy; xref2 = xref;
            %%%% L_BLT_log ~=  log pdf(y|xref)
            LLR_lo=-abs(y2-xref2).^2/N0;  %% Euclidean distance
            LLR_lo=[LLR_lo;zeros(L.^2-ll,M)];
            LLR_lo= LLR_lo(:);
            LLR_lo=LLR_lo';
            %%% add correction
            v2 = abs( N0+2*vVar_PN.*real(xref2.*conj(y2)) );
            VV=imag(conj(xref2).*y2).^2*2.*(vVar_PN/N0)./v2 -0.5*log(v2);
            VV=[VV;zeros(L.^2-ll,M)];
            VV=VV(:);
            VV=VV';
            LLR_log(n,:)=  LLR_lo + VV;
        case 'SP1'
            y2 = yy; xref2 = xref;
            %%%% L_BLT_log ~=  log pdf(y|xref)
            LLR_lo=-abs(abs(y2)-abs(xref2)).^2/N0;  %% Euclidean distance
            LLR_lo=[LLR_lo;zeros(L.^2-ll,M)];
            LLR_lo= LLR_lo(:);
            LLR_lo=LLR_lo';
            delta = angle(y2./xref);
            v2 = abs( N0+2*vVar_PN.*abs(xref2.*y2) );
            VV=0.5*delta.^2./(vVar_PN/N0)./v2 -0.5*log(v2) - 0.5.*delta.^2./vVar_PN;
            VV=[VV;zeros(L.^2-ll,M)];
            VV=VV(:);
            VV=VV';
            LLR_log(n,:)=  LLR_lo + VV;
        case 'SP'
            y2 = yy; xref2 = xref;
            delta = angle(y2./xref); a=N0./(2*abs(xref.*y2).*vVar_PN);
            p_hat = delta./(1+a);
            
            LLR_lo=-abs(y2-xref2.*exp(1j*p_hat)).^2/N0;  %% Euclidean distance
            LLR_lo=[LLR_lo;zeros(L.^2-ll,M)];
            LLR_lo= LLR_lo(:);
            LLR_lo=LLR_lo';
            
            v2 = abs( N0+2*vVar_PN.*abs(xref2.*y2).*cos(p_hat-delta) );
            VV= -0.5*log( v2) - 0.5.*p_hat.^2./vVar_PN;
            VV=[VV;zeros(L.^2-ll,M)];
            VV=VV(:);
            VV=VV';
            
            LLR_log(n,:)=  LLR_lo + VV;
            
        case 'SPit'
            y2 = yy; xref2 = xref;
            delta = angle(y2./xref); a=N0./(2*abs(xref.*y2).*vVar_PN);
            p_hat = delta./(1+a);
            for kk=1:2 %% Newton method to solve saddlepoint equation
                p_hat = p_hat - (sin(p_hat-delta) + p_hat.*a)./(cos(p_hat-delta) + a);
            end
            
            LLR_lo=-abs(y2-xref2.*exp(1j*p_hat)).^2/N0;  %% Euclidean distance
            LLR_lo=[LLR_lo;zeros(L.^2-ll,M)];
            LLR_lo= LLR_lo(:);
            LLR_lo=LLR_lo';
            
            v2 = abs( N0+2*vVar_PN.*abs(xref2.*y2).*cos(p_hat-delta) );
            VV=- 0.5*log( v2) - 0.5.*p_hat.^2./vVar_PN;
            VV=[VV;zeros(L.^2-ll,M)];
            VV=VV(:);
            VV=VV';
            
            LLR_log(n,:)=  LLR_lo + VV;
        case 'EX'
            %% calculation using Gauss-Hermite numerical integration
            nn=15;
            [p_GH, w_GH] = GaussHermite(nn);
            
            gg_fun= @(p) exp(-abs( yy-xref.*exp( 1i*p*sqrt(vVar_PN*2) ) ).^2/N0)/(pi*sqrt(pi)*N0);
            
             LL = zeros(1,L.^2*M);
            for o=1:nn
                VV=gg_fun(p_GH(o))*w_GH(o);
                VV=[VV;zeros(L.^2-ll,M)];
                VV=VV(:);
                VV=VV';
                LL = LL + VV;
            end
            LLR_log(n,:) = log( LL );
        otherwise
            error('unknown integral simplification');
    end
    
end



for k=1:m
    vect1=[];
    vect0=[];
     for f=1:size(outbit1,1)
           vect1=[vect1 [(outbit1(f,k)-1)*L.^2+1:outbit1(f,k)*L.^2]];
           vect0=[vect0 [(outbit0(f,k)-1)*L.^2+1:outbit0(f,k)*L.^2]];
        end
    switch SUM       
        case 'sumexp'
            LLR(:,k)=jac_log_Nclass( LLR_log(:,vect1 ),Coef,Coef_INF,L.^2,M) - jac_log_Nclass( LLR_log(:,vect0 ),Coef,Coef_INF,L.^2,M);
        case 'maxlog'
            LLR(:,k) =( max(LLR_log(:,vect1).')-max(LLR_log(:,vect0).')).';
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