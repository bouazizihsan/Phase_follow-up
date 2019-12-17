function [BCJR]=Phase_track_BCJR_Gauss_App4( y, N0, sigw2, X, pilot_symbols,pilot,Lframe,NClass )


Ns           =length(y);
pilot_pos    = zeros(Ns,1);
pilot_pos(pilot)= pilot_symbols;
M            =length(X);
Len          =NClass-1; % Nombre de class + reste de mixture non classifiées dans unn seul classe
Np           = 2^9;
%Npilots      = length(pilot_symbols);
Phase_limit  = pi; %3*sqrt(sigw2*Lframe/4);  % we suppose the phase varies in the interval (-Phase_limit, Phase_limit)
Threshold    = 3; %pi^2/9;
epsilon      = 3;
dpp          = 2.*Phase_limit./Np;
p           = linspace(-Phase_limit,Phase_limit,Np)';
% ppp=[pp;pp+2.*pi];

pdf_Gauss   = @(t,m,s2) exp( -abs(t-m ).^2./(2*s2) )./sqrt(2*pi*s2); % pdf(t|m), normpdf(x,mu,sigma.^2)
g = fftshift(pdf_Gauss( p ,0 ,sigw2)); %repmat(-dphase_shift(n-1),Np,1 )

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


z1{1,1}=0;
phi1{1,1}=1;
Z_alpha=zeros(Ns,NClass);
Fact_alpha=zeros(Ns,NClass);
Fact_alpha(1)=1;

LT=1; % alpha  class size

%ALPHA

L=1;
alpha_p0=[ 1 Inf 1/(2*pi) ] ;% (moy,var,phi)' ==> (L,3)
parameters1=zeros(Ns,3*NClass);
parameters1(1,1:NClass:2*NClass+1)=[ 1 Inf 1/(2*pi) ];% (moy,var,phi)' ==> (1,3*L)
L_alpha=zeros(Ns,1);
L_alpha(1)=1;
for n=2:Ns
    
    %___________________________________________________
    
% %     if(pilot_pos(n-1)~=0)
% %         
% %         X_Hard=pilot_pos(n-1);
% %         MM=length(X_Hard);
% %         z= z1{1,n-1} + repmat((2/N0)*y(n-1)*conj(X_Hard),LT,1);
% %         
% %         lambda=  exp(abs(z))./sqrt(abs(z)).* exp( -(abs(X_Hard).^2 + abs(y(n-1)).^2 )/N0); %.*exp(real(z)*2/N0);
% %         lambda=phi1{1,n-1}.*lambda;
% %     else
% %         [cel{1,n},uniq]=unique(ii(:,n));
% %         X_Hard=X; %X(ii(uniq,n));
% %         MM=length(X_Hard);
% %         
% %         z=repmat(z1{1,n-1},1,MM) + repmat((2/N0)*y(n-1)*conj(X_Hard),LT,1);
% %         
% %         lambda= exp(abs(z)-repmat(abs(z1{1,n-1}),1,MM)) .* sqrt(repmat(abs(z1{1,n-1}),1,MM)./abs(z)).* repmat(exp(-(abs(X_Hard).^2 + abs(y(n-1)).^2 )/N0),LT,1); % .* exp(real(z)*2/N0);
% %         lambda=repmat(phi1{1,n-1},1,MM).*lambda;
% %     end
% %     
% %     z_gamma =gamma_fct(sigw2,z);
% %     lambda=lambda(:)./sum(lambda(:));
% % 
% %     fig=repmat(lambda./(2*pi*besseli(0,abs(z(:) ))) ,1,Np).*exp(real( repmat(z(:),1,Np).* repmat(exp(-1j*p'),LT*MM,1) ));
% %     fig2=((lambda./(2*pi*besseli(0,abs(z(:) )))).' * exp(real( repmat(z(:),1,Np).* repmat(exp(-1j*p'),LT*MM,1) ))).';
% %     fig=fig';
% %     fig=fig./(repmat(sum(fig),Np,1));
% %     fig2=fig2./(sum(fig2));
% %      figure(1);
% %     plot(p,fig )
% %     hold on
% %     plot(p,fig2, 'r-o')
% %     
% %     hold off
% %     [fact,zT]=KLCluster_U(lambda,z_gamma,epsilon,NClass-1,Threshold);
% %     %[fact,zT]=KLCluster_U2(lambda,z_gamma,NClass,Threshold);
% %     
% %     z1{1,n}=zT;
% %     phi1{1,n}=fact;
% %     alpha = [fact zT];
% %     LT = size(alpha,1);
% % 
% %     Z_alpha(n,1:LT)=zT;
% %     Fact_alpha(n,1:LT)=fact;
% %     VT=1./(abs(zT));    
    %___________________________________________________
    
    if( pilot_pos(n-1)~=0)
        X_Hard=pilot_pos(n-1);
        Mnl_p0{1,n-1} =angle(y(n-1)./X_Hard);
        Snl{1,n-1} = N0./( 2*abs(X_Hard .* y(n-1)) ); % variance
        Snl_Xci{1,n-1} = N0./(2*abs(X_Hard).^2);
        Xci_p0{1,n-1}= pdf_Gauss( abs( y(n-1)./X_Hard) ,1, Snl_Xci{1,n-1} )./ (abs(X_Hard).*sqrt(abs(X_Hard*y(n-1)))) ;
        % %         gamma_p0(:,n)= Xci_p0{1,n}.*pdf_Gauss( p ,Mnl_p0{1,n}, Snl{1,n});
        gamma_exact(:,n-1) =  pdf_Gauss( X_Hard*exp(1j*p) ,y(n-1), N0/2);
    else
        
        Mnl_p0{1,n-1} =angle(y(n-1)./X_Hard);%%   Mnl_p0 = p0 ;
        Snl{1,n-1} = N0./(2*abs(X_Hard).*abs(y(n-1)));
        Snl_Xci{1,n-1} = N0./(2*abs(X_Hard).^2); % variance
        Xci_p0 {1,n-1}= pdf_Gauss( abs( y(n-1)./X_Hard ) ,1, Snl_Xci{1,n-1})./ (sqrt(abs(X_Hard*y(n-1))).*abs(X_Hard) );
    end
    
    
    PHI_nl_len_p0=zeros(L,length(Xci_p0{1,n-1}) );
    Snl_tilde =zeros(L,length(Xci_p0{1,n-1}) );
    Mnl_tilde_p0=zeros(L,length(Xci_p0{1,n-1}) );
    for k=1:L
        if(alpha_p0(k+L)==Inf)
            PHI_nl_len_p0(k,:)=alpha_p0(2*L+k).* Xci_p0{1,n-1}./(2*pi);
            Mnl_tilde_p0(k,:)= Mnl_p0{1,n-1};
            Snl_tilde(k,:)=  Snl{1,n-1};
        else
            PHI_nl_len_p0(k,:)=alpha_p0(2*L+k).* Xci_p0{1,n-1} .* pdf_Gauss( alpha_p0(k) ,Mnl_p0{1,n-1},alpha_p0(k+L)+ Snl{1,n-1});
            Mnl_tilde_p0(k,:)=(alpha_p0(k).*Snl{1,n-1} + Mnl_p0{1,n-1}.*alpha_p0(L+k))./(alpha_p0(L+k) + Snl{1,n-1});
            Snl_tilde(k,:)= alpha_p0(L+k) .* Snl{1,n-1} ./ (alpha_p0(L+k) + Snl{1,n-1});
            
        end
    end
    
    
    %     end
    PHI_nl_len_p0=PHI_nl_len_p0/sum(PHI_nl_len_p0(:));
    %[PHI_p0_sort,kk_sort]=sort(PHI_p0);
    a=PHI_nl_len_p0.';
   
    PHI_nl_len_p=a(:);
    %Len=length(PHI_nl_len_p);
    b=Mnl_tilde_p0.';

    Mnl_tilde_p=b(:);
    
    c=Snl_tilde.' +sigw2;
    Snl_tild=c(:);
    
% %     fig=repmat(PHI_nl_len_p.',Np,1).*pdf_Gauss( repmat(p,1,L*length(Xci_p0{1,n-1})) ,repmat(Mnl_tilde_p.',Np,1), repmat(Snl_tild.',Np,1));
% %     fig2=(PHI_nl_len_p.')*(pdf_Gauss( repmat(p,1,L*length(Xci_p0{1,n-1})) ,repmat(Mnl_tilde_p.',Np,1), repmat(Snl_tild.',Np,1))).';
% %     figure(2);
% %     
% %     plot(p,fig./(repmat(sum(fig),Np,1)) )
% %     hold on
% %     plot(p,fig2'./(sum(fig2')),'r-o' )
% %     
% %     hold off

    alph=[];
    j=1;
    while ((j <= Len) && (norm(PHI_nl_len_p)>0))
        [Mnl,S]=Gaussian1(Mnl_tilde_p,Snl_tild,PHI_nl_len_p);
        
%         [var, indx] = max(PHI_nl_len_p);
%         f=DKL(Mnl_tilde_p(indx),Mnl_tilde_p,Snl_tild(indx),Snl_tild);
        f=DKL(Mnl,Mnl_tilde_p,S,Snl_tild);
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
        [Mnl,SR]=Gaussian1(Mnl_tilde_p,Snl_tild,PHI_nl_len_p);
        betaj= sum(PHI_nl_len_p);
        if(SR< Threshold)
            alph= [alph [Mnl ;SR;betaj]];
        else
            alph= [alph [Mnl ;Inf;betaj]];
        end
    end
    
    variab=alph.';
    variab=variab(:);
    alpha_p0=variab';
    L=length(alpha_p0)./3;
    L_alpha(n)=L;
    variab=[alph zeros(3,NClass-L)];
    variab=variab.';
    variab=variab(:);
    parameters1(n,:)=variab.';
   for k=1:L 
    Palpha= sqrt(alpha_p0(L+k)).*randn + alpha_p0(k) ;
   end
   YY=y(n).*exp(-1j*Palpha);
   [X_Hard,i]=Hard_Decision(YY,X);
    %% calcul exact
% %     
% %         alpha_exact(:,n) = ifft( fft( alpha_exact(:,n-1).* gamma_exact(:,n-1), Np ) .* fft(g, Np)  );
% %         alpha_p0_Len(:,n)=(alpha_p0(2*L+1:3*L)*pdf_Gauss( repmat(p,1,L) ,repmat(alpha_p0(1:L),Np,1), repmat(alpha_p0(L+1:2*L),Np,1)).')';
% %         alpha_plot(:,n)= exp(real( repmat(z1{1,n},1,Np).* repmat(exp(-1j*p'),LT,1)   ))' * (phi1{1,n}./(2*pi*besseli(0,abs(z1{1,n}) ) ) );% + fact_inf./(2*pi); z=0 --> pdf=1./(2*pi)
% % %         alpha_plotG(:,n)= pdf_Gauss( repmat(p.',LT,1) ,repmat(angle(z1{1,n}),1,Np), repmat(1./abs(z1{1,n}),1,Np))' * phi1{1,n} ;% + fact_inf./(2*pi); z=0 --> pdf=1./(2*pi)
% %        
% %         fig=repmat(alpha_p0(2*L+1:3*L),Np,1).*pdf_Gauss( repmat(p,1,L) ,repmat(alpha_p0(1:L),Np,1), repmat(alpha_p0(L+1:2*L),Np,1));
% %         if(norm(PHI_nl_len_p))
% %         fig2=repmat(PHI_nl_len_p.',Np,1).*pdf_Gauss( repmat(p,1,length(Mnl_tilde_p.')),repmat(Mnl_tilde_p.',Np,1), repmat(Snl_tild.',Np,1));
% %         end
% % %         normalize
% %         alpha_exact(:,n) = alpha_exact(:,n)./sum( alpha_exact(:,n));   
% %         alpha_p0_Len(:,n)=alpha_p0_Len(:,n)./(sum(alpha_p0_Len(:,n)));
% %         alpha_plot(:,n)=alpha_plot(:,n)./sum(alpha_plot(:,n));
% % %         alpha_plotG(:,n)=alpha_plotG(:,n)./sum(alpha_plotG(:,n));
% %     
% %         figure(3)
% %         plot(p,alpha_p0_Len(:,n),'b')
% %         hold on
% %         plot(p,alpha_exact(:,n), 'g')
% %         plot(p,alpha_plot(:,n),'r') 
% % %         plot(p,alpha_plotG(:,n),'r')
% %         plot(p,fig./(repmat(sum(fig),Np,1)))
% %         if(norm(PHI_nl_len_p))
% %          hold on
% %          plot(p,fig2./(repmat(sum(fig2),Np,1)),'--' )
% %          end
% %         title({'alpha, ',n})
% %         hold off
% %     
% %     
% %     moy=(p')*( alpha_plot(:,n)./(sum(alpha_plot(:,n))));
% %     var=(p').^2*(alpha_plot(:,n)./(sum(alpha_plot(:,n),1)))- moy.^2;
% %     
% % %     moyG=(p')*( alpha_plotG(:,n)./(sum(alpha_plotG(:,n))));
% % %     varG=(p').^2*(alpha_plotG(:,n)./(sum(alpha_plotG(:,n),1)))- moyG.^2;
% % %     
% % %     (alpha_p0(L+1:2*L))
% % %     (VT')
% % %     (var)
% % %     varG
end


%beta
L=1;
beta_p0=[1 Inf 1/(2*pi)];
parameters2=zeros(Ns,3*NClass);
parameters2(Ns,1:NClass:2*NClass+1)=[ 1 Inf 1/(2*pi)];
L_beta=zeros(Ns,1);
L_beta(Ns)=1;
for n=Ns-1:-1:1
    
    if( pilot_pos(n+1)~=0)
        X_Hard=pilot_pos(n+1);
        Mnl_p0{1,n+1} =angle(y(n+1)./X_Hard);
        Snl{1,n+1} = N0./( 2*abs(X_Hard .* y(n+1)) ); % variance
        Snl_Xci{1,n+1} = N0./(2*abs(X_Hard).^2);
        Xci_p0{1,n+1}= pdf_Gauss( abs( y(n+1)./X_Hard) ,1, Snl_Xci{1,n+1} )./ (abs(X_Hard).*sqrt(abs(X_Hard*y(n+1)))) ;
        % %         gamma_p0(:,n)= Xci_p0{1,n}.*pdf_Gauss( p ,Mnl_p0{1,n}, Snl{1,n});
        gamma_exact(:,n+1) =  pdf_Gauss( X_Hard*exp(1j*p) ,y(n+1), N0/2);
    else
        
        Mnl_p0{1,n+1} =angle(y(n+1)./X_Hard);%%   Mnl_p0 = p0 ;
        Snl{1,n+1} = N0./(2*abs(X_Hard).*abs(y(n+1)));
        Snl_Xci{1,n+1} = N0./(2*abs(X_Hard).^2); % variance
        Xci_p0 {1,n+1}= pdf_Gauss( abs( y(n+1)./X_Hard ) ,1, Snl_Xci{1,n+1})./ (sqrt(abs(X_Hard*y(n+1))).*abs(X_Hard) );
    end
    
    
    PHI_nl_len_p0=zeros(L,length(Xci_p0{1,n+1}) );
    Snl_tilde =zeros(L,length(Xci_p0{1,n+1}) );
    Mnl_tilde_p0=zeros(L,length(Xci_p0{1,n+1}) );
    
    for k=1:L
        
        if(beta_p0(k+L)==Inf)
            PHI_nl_len_p0(k,:)=beta_p0(2*L+k).* Xci_p0{1,n+1}./(2*pi);
            Mnl_tilde_p0(k,:)= Mnl_p0{1,n+1};
            Snl_tilde(k,:)=  Snl{1,n+1};
        else
            PHI_nl_len_p0(k,:)=beta_p0(2*L+k).* Xci_p0{1,n+1} .* pdf_Gauss( beta_p0(k) ,Mnl_p0{1,n+1},beta_p0(k+L)+ Snl{1,n+1});
            Mnl_tilde_p0(k,:)=(beta_p0(k).*Snl{1,n+1} + Mnl_p0{1,n+1}.*beta_p0(L+k))./(beta_p0(L+k) + Snl{1,n+1});
            Snl_tilde(k,:)= beta_p0(L+k) .* Snl{1,n+1} ./ (beta_p0(L+k) + Snl{1,n+1});
            
        end
        
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
    while ((j <= Len) && (norm(PHI_nl_len_p)>0))
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
        if(S<Threshold)
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
    variab=[beta zeros(3,NClass-L)];
    variab=variab.';
    variab=variab(:);
    parameters2(n,:)=variab;
    
    for k=1:L 
    Pbeta= sqrt(beta_p0(L+k)).*randn + beta_p0(k) ;
   end
   YY=y(n).*exp(-1j*Pbeta);
   [X_Hard,i]=Hard_Decision(YY,X);
    
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

Pmean = ( p'*p_apo )*dpp;
Pvar  = ( (p.^2)'*p_apo )*dpp - Pmean.^2;
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

dist=abs(B-A);
[ind v]=find(dist>pi);
dist(ind)=2*pi-dist(ind);
f=dist.^2./C + D./C +0.5*log(C./D) - 1;

end

function z=gamma_fct(sig,z)

z=z./(1+abs(z).*sig);

end