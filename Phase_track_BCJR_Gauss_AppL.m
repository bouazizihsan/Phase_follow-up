function [BCJR]=Phase_track_BCJR_Gauss_AppL( y, N0, sigw2, X, pilot_symbols,pilot,Lframe,NClass )

epsilon=3;
Ns           =length(y);
pilot_pos    = zeros(Ns,1);
pilot_pos(pilot)= pilot_symbols;
M            =length(X);
Len            = NClass-1;
Np           = 2^9;
%Npilots      = length(pilot_symbols);
Phase_limit  = 3*sqrt(sigw2*Lframe/4);  % we suppose the phase varies in the interval (-Phase_limit, Phase_limit)
dpp          = 2.*Phase_limit./Np;
pp           = linspace(-Phase_limit,Phase_limit,Np)';
% ppp=[pp;pp+2.*pi];

pdf_Gauss   = @(t,m,s2) exp( -abs(t-m ).^2./(2*s2) )./sqrt(2*pi*s2); % pdf(t|m), normpdf(x,mu,sigma.^2)
g = fftshift(pdf_Gauss( pp ,0 ,sigw2)); %repmat(-dphase_shift(n-1),Np,1 )

%processing
phase=angle(y(pilot)./pilot_symbols);
dephasage= diff(phase);
fi=find(dephasage > pi); dephasage(fi)= dephasage(fi)-2.*pi;
fi=find(dephasage < -pi); dephasage(fi)= dephasage(fi)+2.*pi;
new_phases= phase(1) +[0;cumsum(dephasage)];
phase_shift=interp1(pilot', new_phases', (1:Ns)');
dphase_shift=diff(phase_shift);


y=y.*exp(-1j.*phase_shift);

gamma_exact = zeros(Np,Ns);
gamma_p0 = zeros(Np,Ns);

alpha_exact=zeros(Np,Ns);
beta_exact=zeros(Np,Ns);

alpha_p0_Len =  zeros(Np,Ns);
beta_p0_Len =  zeros(Np,Ns);

alpha_p0_Len(:,1)=1;
alpha_exact(:,1)=1;
beta_exact(:,Ns)=1;
beta_p0_Len(:,Ns)=1;

Xci_p0=zeros(Ns,M);
Mnl_p0=zeros(Ns,M);
Snl=zeros(Ns,M);
Snl_Xci=zeros(Ns,M);
for n=1:Ns
    if( pilot_pos(n)~=0)       
        Mnl_p0(n,1) =angle(y(n)./pilot_pos(n));        
        Snl(n,1) = N0./(2*abs(pilot_pos(n)).*abs(y(n)) ); % variance
        Snl(n,2:M)=Inf;
        Snl_Xci(n,1) = N0./(2*abs(pilot_pos(n)).^2);
        Xci_p0(n,1)= pdf_Gauss( abs( y(n)./pilot_pos(n)) ,1, Snl_Xci(n,1) )./ (abs(pilot_pos(n)).*sqrt(abs(pilot_pos(n)*y(n)))) ;
        gamma_p0(:,n)= Xci_p0(n,1).*pdf_Gauss( pp ,Mnl_p0(n,1), Snl(n,1)); 
        gamma_exact(:,n) =  pdf_Gauss( pilot_pos(n)*exp(1j*pp) ,y(n), N0/2);        
    else
        for nn=1:M
            Mnl_p0(n,nn) =angle(y(n)./X(nn));%%   Mnl_p0 = p0 ; 
            Snl(n,nn) = N0./(2*abs(X(nn)).*abs(y(n)));
            Snl_Xci(n,nn) = N0./(2*abs(X(nn)).^2); % variance
            %Xci (n,nn)= pdf_Gauss( real( y(n)./X(nn) ) ,1, Snl(n,nn) )./ abs(X(nn)).^2 ;
            Xci_p0 (n,nn)= pdf_Gauss( abs( y(n)./X(nn) ) ,1, Snl_Xci(n,nn) )./ (sqrt(abs(X(nn)*y(n))).*abs(X(nn)) );
            gamma_p0(:,n) = gamma_p0(:,n) + Xci_p0(n,nn).*pdf_Gauss( pp ,Mnl_p0(n,nn), Snl(n,nn));                    
            gamma_exact(:,n) = gamma_exact(:,n) + pdf_Gauss( X(nn)*exp(1j*pp) ,y(n), N0/2);
            
        end
        
    end
end

%ALPHA

L=1;
alpha_p0=[ 1 Inf 1/(2*pi) ] ;% (moy,var,phi)' ==> (L,3)
fact=zeros(1,Ns);
fact(1)=1;
parameters1=zeros(Ns,3*NClass);
parameters1(1,1:NClass:2*NClass+1)=[ 1 Inf 1/(2*pi) ];% (moy,var,phi)' ==> (1,3*L)
L_alpha=zeros(Ns,1);
L_alpha(1)=1;
for n=2:Ns
    PHI_nl_len_p0=zeros(L,M);
    Snl_tilde =zeros(L,M);
    Mnl_tilde_p0=zeros(L,M);
 
    for k=1:L
        
        if(alpha_p0(k+L)==Inf)
            PHI_nl_len_p0(k,:)=alpha_p0(2*L+k).* Xci_p0(n-1,:)./(2*pi);
            Mnl_tilde_p0(k,:)= Mnl_p0(n-1,:);
            Snl_tilde(k,:)=  Snl(n-1,:);
        else
            PHI_nl_len_p0(k,:)=alpha_p0(2*L+k).* Xci_p0(n-1,:) .* pdf_Gauss( alpha_p0(k) ,Mnl_p0(n-1,:),alpha_p0(k+L)+ Snl(n-1,:));
            Mnl_tilde_p0(k,:)=(alpha_p0(k).*Snl(n-1,:) + Mnl_p0(n-1,:).*alpha_p0(L+k))./(alpha_p0(L+k) + Snl(n-1,:));
            Snl_tilde(k,:)= alpha_p0(L+k) .* Snl(n-1,:) ./ (alpha_p0(L+k) + Snl(n-1,:));
            
        end
        
    end
    if (pilot_pos(n-1)~=0)
        PHI_nl_len_p0(:,2:end)=[];
        Mnl_tilde_p0(:,2:end)=[];
        Snl_tilde(:,2:end)=[];
    end
    %     end
    PHI_nl_len_p0=PHI_nl_len_p0/sum(PHI_nl_len_p0(:));   
    %[PHI_p0_sort,kk_sort]=sort(PHI_p0);
    a=PHI_nl_len_p0.';
    aa=a(:);
    PHI_nl_len_p=a(:);
    %Len=length(PHI_nl_len_p);
    b=Mnl_tilde_p0.';
    bb=b(:);
    Mnl_tilde_p=b(:);
    
    c=Snl_tilde.' +sigw2;
    cc=c(:);
    Snl_tild=c(:);
    
   
    alph=[];
    j=1;
    while ((j <= NClass-1) && (norm(PHI_nl_len_p)>0))
        [var, indx] = max(PHI_nl_len_p);
        
        f=DKL(Mnl_tilde_p(indx),Mnl_tilde_p,Snl_tild(indx),Snl_tild);
        aux=find(f<epsilon);
        sum_phi=sum(PHI_nl_len_p(aux));
        %la moyenne
        [Mnl,S]=Gaussian1(Mnl_tilde_p(aux),Snl_tild(aux),PHI_nl_len_p(aux));
        alph= [alph [Mnl ;S;sum_phi]];

        Mnl_tilde_p(aux)= [];
        Snl_tild(aux)=[];
        PHI_nl_len_p(aux)=[];
         j=j+1;
    end
    if(norm(PHI_nl_len_p))
        [Mnl,S]=Gaussian1(Mnl_tilde_p,Snl_tild,PHI_nl_len_p);
        betaj= sum(PHI_nl_len_p);
        if(S<3)
            alph= [alph [Mnl ;S;betaj]];
        else
            alph= [alph [Mnl ;Inf;betaj]];
        end
    end
        
    variab=alph.';
    variab=variab(:);
    alpha_p0=variab';
    L=length(alpha_p0)./3;
    L_alpha(n)=L;
    alp=[alph zeros(3,Len+1-L)];
    alp=alp.';
    alp=alp(:);
    parameters1(n,:)=alp.';
    
    
    %% calcul exact
    
% %     alpha_exact(:,n) = ifft( fft( alpha_exact(:,n-1).* gamma_exact(:,n-1), Np ) .* fft(g, Np)  );
% %     alpha_p0_Len(:,n)=(alpha_p0(2*L+1:3*L)*pdf_Gauss( repmat(pp,1,L) ,repmat(alpha_p0(1:L),Np,1), repmat(alpha_p0(L+1:2*L),Np,1)).')';
% %     
% %     %normalize
% %     alpha_exact(:,n) = alpha_exact(:,n)./sum( alpha_exact(:,n));
% %     
% %     alpha_p0_Len(:,n)=alpha_p0_Len(:,n)./(sum(alpha_p0_Len(:,n)));
% %     
% %     
% %     figure(2)
% %     plot(pp,alpha_p0_Len(:,n),'b')
% %     hold on
% %     
% %     plot(pp,alpha_exact(:,n), 'g')
% %     hold on
% %     som=aa'*pdf_Gauss( repmat(pp,1,length(bb)),repmat(bb',Np,1),repmat(cc',Np,1)).';
% %     som=som./sum(som);
% %     plot(pp,som,'k')
% % 
% %     title({'alpha, ',n})
% %     hold off
% %     
% %     
    % moy=(pp')*( alpha_item./repmat(sum(alpha_item,1),Np,1));
    % var=(pp').^2*(alpha_item./repmat(sum(alpha_item,1),Np,1))- moy.^2;
    
    
end


%beta
L=1;
beta_p0=[1 Inf 1/(2*pi)];
parameters2=zeros(Ns,3*NClass);
parameters2(Ns,1:NClass:2*NClass+1)=[ 1 Inf 1/(2*pi)];
L_beta=zeros(Ns,1);
L_beta(Ns)=1;
for n=Ns-1:-1:1
    
    PHI_nl_len_p0=zeros(L,M);
    Snl_tilde =zeros(L,M);
    Mnl_tilde_p0=zeros(L,M);
 
    for k=1:L
        
        if(beta_p0(k+L)==Inf)
            PHI_nl_len_p0(k,:)=beta_p0(2*L+k)/(2*pi).* Xci_p0(n+1,:);
            Mnl_tilde_p0(k,:)= Mnl_p0(n+1,:);
            Snl_tilde(k,:)=  Snl(n+1,:);
        else
            PHI_nl_len_p0(k,:)=beta_p0(2*L+k).* Xci_p0(n+1,:) .* pdf_Gauss( beta_p0(k) ,Mnl_p0(n+1,:),beta_p0(k+L)+ Snl(n+1,:));
            Mnl_tilde_p0(k,:)=(beta_p0(k).*Snl(n+1,:) + Mnl_p0(n+1,:).*beta_p0(L+k))./(beta_p0(L+k) + Snl(n+1,:));
            Snl_tilde(k,:)= beta_p0(L+k) .* Snl(n+1,:) ./ (beta_p0(L+k) + Snl(n+1,:));
            
        end
        
    end
    if (pilot_pos(n+1)~=0)
        PHI_nl_len_p0(:,2:end)=[];
        Mnl_tilde_p0(:,2:end)=[];
        Snl_tilde(:,2:end)=[];
    end
    %     end
    PHI_nl_len_p0=PHI_nl_len_p0/sum(PHI_nl_len_p0(:));
    
    %[PHI_p0_sort,kk_sort]=sort(PHI_p0);
    a=PHI_nl_len_p0.';
    aa=a(:);
    PHI_nl_len_p=a(:);
    %Len=length(PHI_nl_len_p);
    b=Mnl_tilde_p0.';
    bb=b(:);
    Mnl_tilde_p=b(:);
    
    c=Snl_tilde.' +sigw2;
    cc=c(:);
    Snl_tild=c(:);
    
%     figure(3)
%     plot(pp,repmat(PHI_nl_len_p,1,Np).*pdf_Gauss( repmat(pp,1,length(Mnl_tilde_p)),repmat(Mnl_tilde_p',Np,1),repmat(Snl_tild',Np,1)).')
    
    betaj=[];
    beta=[];
    j=1;
    while ((j <= NClass-1) && (norm(PHI_nl_len_p)>0))
        [var, indx] = max(PHI_nl_len_p);
        
        f=DKL(Mnl_tilde_p(indx),Mnl_tilde_p,Snl_tild(indx),Snl_tild);
        aux=find(f<epsilon);
        sum_phi=sum(PHI_nl_len_p(aux));
        betaj=[betaj sum_phi];
        %la moyenne
        [Mnl,S]=Gaussian1(Mnl_tilde_p(aux),Snl_tild(aux),PHI_nl_len_p(aux));
        beta= [beta [Mnl ;S;sum_phi]];

        Mnl_tilde_p(aux)= [];
        Snl_tild(aux)=[];
        PHI_nl_len_p(aux)=[];
        
        j=j+1;
    end
    
    if(norm(PHI_nl_len_p))
        [Mnl,S]=Gaussian1(Mnl_tilde_p,Snl_tild,PHI_nl_len_p);
        betaj= sum(PHI_nl_len_p);
        if(S<3)
            beta= [beta [Mnl ;S;betaj]];
        else
            beta= [beta [Mnl ;Inf;betaj]];
        end
              
    end
    variab=beta';
    variab=variab(:);
    beta_p0=variab';
    L=length(beta_p0)./3;
    L_beta(n)=L;
    bet=[beta zeros(3,NClass-L)];
    bet=bet.';
    bet=bet(:);
    parameters2(n,:)=bet;
    
    
% %     %% calcul exact
% %     
% %     beta_exact(:,n) = ifft( fft( beta_exact(:,n+1).* gamma_exact(:,n+1), Np ) .* fft(g, Np)  );
% %     beta_p0_Len(:,n)=beta_p0(2*L+1:3*L)*pdf_Gauss( repmat(pp,1,L) ,repmat(beta_p0(1:L),Np,1), repmat(beta_p0(L+1:2*L),Np,1)).';
% %     
% %     %normalize
% %     beta_exact(:,n) = beta_exact(:,n)./sum( beta_exact(:,n));
% %     
% %     beta_p0_Len(:,n)=beta_p0_Len(:,n)./(sum(beta_p0_Len(:,n)));
% %     
% %     
% %     figure(2)
% %     plot(pp,beta_p0_Len(:,n),'b')
% %     hold on
% %     
% %     plot(pp,beta_exact(:,n), 'g')
% %     hold on
% %     som=aa'*pdf_Gauss( repmat(pp,1,length(bb)),repmat(bb',Np,1),repmat(cc',Np,1)).';
% %     som=som./sum(som);
% %     plot(pp,som,'k')
% % 
% %     title({'BETA--> ',n})
% %     hold off
% %     
       
end


p_apo=beta_exact.*alpha_exact;
p_apo=p_apo./repmat(sum(p_apo,1),Np,1)/dpp;

Pmean = ( pp'*p_apo )*dpp;
Pvar  = ( (pp.^2)'*p_apo )*dpp - Pmean.^2;
Pmean=Pmean.';
Pvar=Pvar.';

% Pmean(pilot)=[];
% Pvar(pilot)=[];
parameters1(pilot,:)=[];
parameters2(pilot,:)=[];
L_alpha(pilot,:)=[];
L_beta(pilot,:)=[];

BCJR.parameters1=parameters1;
BCJR.parameters2=parameters2;
BCJR.L_alpha=L_alpha;
BCJR.L_beta=L_beta;
BCJR.Pmean=Pmean;
BCJR.Pvar=Pvar;
BCJR.phase_shift=phase_shift;
end

function f=DKL(A,B,C,D)
%f= sum(A.*log(A./repmat(B,1,size(A,2))),1);
dist=abs(B-A);
[ind v]=find(dist>pi);
dist(ind)=2*pi-dist(ind);

f=dist.^2./C + D./C +0.5*log(C./D) - 1;

end