function [BCJR]=Phase_track_Tikhonov_Hard( y, N0, sigw2, X, pilot_symbols,pilot,Lframe,NClass )

epsilon      =12;
Ns           =length(y);
pilot_pos(pilot)= pilot_symbols;
M            =length(X);
Np           = 2^9;
% BCJR.alpha=zeros(Np,Ns);
Npilots      = length(pilot_symbols);
Phase_limit  = pi; %3*sqrt(sigw2*Lframe/4);  % we suppose the phase varies in the interval (-Phase_limit, Phase_limit)
dpp          = 2.*Phase_limit./Np;
pp           = linspace(-Phase_limit,Phase_limit,Np)';
Threshold    = Phase_limit^2/9;

pdf_Gauss   = @(t,m,s2) exp( -abs(t-m ).^2./(2*s2) )./sqrt(2*pi*s2); % pdf(t|m), normpdf(x,mu,sigma.^2)
g =fftshift(pdf_Gauss( pp ,0 ,sigw2));

%processing
phase=angle(y(pilot)./pilot_symbols);
dephasage= diff(phase);
fi=find(dephasage> pi); dephasage(fi)= dephasage(fi)-2.*pi;
fi=find(dephasage< -pi); dephasage(fi)= dephasage(fi)+2.*pi;
new_phases= phase(1) +[0;cumsum(dephasage)];
phase_shift=interp1( pilot' , new_phases, 1:Ns);
dphase_shift=diff(phase_shift);

y=y.*exp(-1j.*phase_shift');
yy=repmat(y.',Np,1);
p=repmat(pp,1,Ns);
YY=yy.*exp(-1j.*p);

[XX,i]=Hard_Decision(YY(:),X); %% Hard_Decision(YY(:),X); hard_QAM
XX=reshape(XX,Np,Ns);
ii=reshape(i,Np,Ns);

alpha_plot=zeros(Np,Ns);
alpha_exact=zeros(Np,Ns);
alpha_exact(:,1)=1/(2*pi);
alpha_exact(:,1)=alpha_exact(:,1)./sum(alpha_exact(:,1));

beta_plot=zeros(Np,Ns);
beta_exact=zeros(Np,Ns);
beta_exact(:,Ns)=1/(2*pi);
beta_exact(:,Ns)=beta_exact(:,Ns)./sum(beta_exact(:,Ns));

gamma=zeros(Np,Ns);

% % for n=1:Ns
% %     if(pilot_pos(n)~=0)
% %         gamma(:,n) = pdf_Gauss( y(n) ,pilot_pos(n).*exp(1j*pp), N0/2);
% %     else
% %         for nn=1:M
% %             gamma(:,n) = gamma(:,n) + pdf_Gauss( y(n) , x(nn).*exp(1j*pp), N0/2);
% %         end
% %     end
% %     gamma(:,n)=gamma(:,n)./sum(gamma(:,n));
% % end

z1{1,1}=0;
phi1{1,1}=1;
Z_alpha=zeros(Ns,NClass);
Fact_alpha=zeros(Ns,NClass);
Fact_alpha(1)=1;

L_alpha=zeros(Ns,1);
L_alpha(1)=1;

L=1; % alpha  class size
% Ll=1;

for n=2:Ns
    
    if(pilot_pos(n-1)~=0)
        
        X_Hard=pilot_pos(n-1);
    else
        [cel{1,n},uniq]=unique(ii(:,n));
        X_Hard=X; %X(ii(uniq,n));
    end
    MM=length(X_Hard);
    z=repmat(z1{1,n-1},1,MM) + repmat((2/N0)*y(n-1)*conj(X_Hard),L,1);
    
    k=abs(repmat(z1{1,n-1},1,MM) ) + abs(repmat((2/N0)*y(n-1)*conj(X_Hard),L,1));
    var=repmat(abs(z1{1,n-1}),1,MM);
    [l,c,v]=find(var==0);
    
    u=(angle(repmat(z1{1,n-1},1,MM))./abs(repmat((2/N0)*y(n-1)*conj(X_Hard),L,1)) + angle(repmat((2/N0)*y(n-1)*conj(X_Hard),L,1))./ abs(repmat(z1{1,n-1},1,MM)) ) ./( 1./abs(repmat(z1{1,n-1},1,MM)) + 1./abs(repmat((2/N0)*y(n-1)*conj(X_Hard),L,1))) ;
    
    lambda= log(repmat(phi1{1,n-1},1,MM))+(abs(z)-repmat(abs(z1{1,n-1}),1,MM))-0.5*log(abs(z)./repmat(abs(z1{1,n-1}),1,MM))-repmat(abs(X_Hard).^2/N0,L,1);
    
    if(l)
        u(l,c)=(angle(var(l,c))./abs(repmat((2/N0)*y(n-1)*conj(X_Hard),L,1)) ) ./(1./abs(repmat((2/N0)*y(n-1)*conj(X_Hard),L,1))) ;
        var=repmat(phi1{1,n-1},1,MM);
        lambda(l,c)=log(var(l,c))+abs(z(l,c))-0.5*log(abs(z(l,c)))-(abs(X_Hard).^2)/N0;
    end
    lambda=lambda(:);
    lambda= exp(lambda-max(lambda)-log( sum( exp(lambda-max(lambda))) )  );
   
%     z_gamma=k.*exp(1j*u);
    z_gamma =gamma_fct(sigw2,z(:));
    [fact,zT]=KLCluster_U(lambda,z_gamma,epsilon,NClass-1,Threshold);
    %     [fact,zT]=KLCluster_U2(lambda,z_gamma,NClass,Threshold);
    
    z1{1,n}=zT;
    phi1{1,n}=fact;
    L = length(zT);   
    Z_alpha(n,1:L)=zT;
    Fact_alpha(n,1:L)=fact;
    L_alpha(n)=L;
    
    
    %%%_____________ plot Alpha __________
    
    % %         alpha_plot(:,n)= exp(real( repmat(alpha(:,2),1,Np).* repmat(exp(-1j*pp'),L,1)   ))' * (alpha(:,1)./(2*pi*besseli(0,abs(alpha(:,2)) ) ) );% + fact_inf./(2*pi); z=0 --> pdf=1./(2*pi)
    % %         alpha_plot(:,n)=alpha_plot(:,n)./sum(alpha_plot(:,n));
    % %
    % %         alpha_exact(:,n)=ifft(fft(alpha_exact(:,n-1).*gamma(:,n-1)).*fft(g));
    % %         alpha_exact(:,n)=alpha_exact(:,n)./sum(alpha_exact(:,n));
    % %
    % %
    % %         alpha_tkh(:,n)= exp(real( repmat(z_G,1,Np).* repmat(exp(-1j*pp'),Ll,1)   ))' * (factG./(2*pi*besseli(0,abs(z_G)) ) );% + fact_inf./(2*pi); z=0 --> pdf=1./(2*pi)
    % %         alpha_tkh(:,n)=alpha_tkh(:,n)./sum(alpha_tkh(:,n));
    % %
    % %
    % %
    % %         plot(pp,alpha_plot(:,n), 'r');
    % %         hold on;
    % %         plot(pp,alpha_exact(:,n), 'g');
    % %
    % %         plot(pp,alpha_tkh(:,n), 'k');
    % %         xlabel({'n= ',n});
    % %         hold off
    % %     legend('Tikh','exact','Gauss','Tikh-G')
    %%%_____________ END plot Alpha __________
end

%Beta
beta=[1 0];
z2{1,Ns}=0;
phi2{1,Ns}=1;

L=1;
Z_beta=zeros(Ns,NClass);
Fact_beta=zeros(Ns,NClass);
Fact_beta(Ns)=1;

L_beta=zeros(Ns,1);
L_beta(Ns)=1;

for n=Ns-1:-1:1
    if(pilot_pos(n+1) ~= 0)
        X_Hard=pilot_pos(n+1);
    else
        [cel{1,n},uniq]=unique(ii(:,n));
        X_Hard=X; %X(ii(uniq,n)); %% X_Hard = X(cel{1,n})
    end
    
    MM=length(X_Hard);
    z=repmat(z2{1,n+1},1,MM) + repmat((2/N0)*y(n+1)*conj(X_Hard),L,1);
    k=abs(repmat(z2{1,n+1},1,MM)) + abs(repmat((2/N0)*y(n+1)*conj(X_Hard),L,1));
    var=repmat(abs(z2{1,n+1}),1,MM);
    [l,c,v]=find(var==0);
    u=(angle(repmat(z2{1,n+1},1,MM))./(repmat((2/N0)*y(n+1)*conj(X_Hard),L,1)) + angle(repmat((2/N0)*y(n+1)*conj(X_Hard),L,1))./abs(repmat(z2{1,n+1},1,MM)) ) ./( 1./abs(repmat(z2{1,n+1},1,MM)) + 1./abs(repmat((2/N0)*y(n+1)*conj(X_Hard),L,1))) ;
    lambda= log(repmat(phi2{1,n+1},1,MM))+(abs(z)-repmat(abs(z2{1,n+1}),1,MM))-0.5*log(abs(z)./repmat(abs(z2{1,n+1}),1,MM))-repmat(abs(X_Hard).^2/N0,L,1);
    
    if(l)
        u(l,c)=(angle(var(l,c))./(repmat((2/N0)*y(n+1)*conj(X_Hard),L,1)) ) ./(1./abs(repmat((2/N0)*y(n+1)*conj(X_Hard),L,1))) ;
        var=repmat(phi2{1,n+1},1,MM);
        lambda(l,c)=log(var(l,c))+abs(z(l,c))-0.5*log(abs(z(l,c)))-(abs(X_Hard).^2)/N0;
    end
    lambda=lambda(:);
    lambda= exp(lambda-max(lambda)-log( sum( exp(lambda-max(lambda))) )  );
    if(sum(lambda)==0)
      st=1;  
    end
%     zg=k.*exp(1j*u);
    z_gamma=gamma_fct(sigw2,z(:));
    [fact,zT]=KLCluster_U(lambda,z_gamma,epsilon,NClass-1,Threshold);
    %     [fact,zT]=KLCluster_U2(lambda,z_gamma,NClass,Threshold);
    
    z2{1,n}=zT;
    phi2{1,n}=fact;
    
    L=length(zT);
    Z_beta(n,1:L)=zT;
    Fact_beta(n,1:L)=fact;
    L_beta(n)=L;
   
    
    %%%_____________ Plot BETA __________
    % %
    % %         beta_plot(:,n)= exp(real( repmat(beta(:,2),1,Np).* repmat(exp(-1j*pp'),L,1)   ))' * (beta(:,1)./(2*pi*besseli(0,abs(beta(:,2)) ) ) );% + fact_inf./(2*pi); z=0 --> pdf=1./(2*pi)
    % %         beta_plot(:,n)=beta_plot(:,n)./sum(beta_plot(:,n));
    % %
    % %         beta_exact(:,n)=ifft(fft(beta_exact(:,n+1).*gamma(:,n+1)).*fft(g));
    % %         beta_exact(:,n)=beta_exact(:,n)./sum(beta_exact(:,n));
    % %
    % %         beta_tkh(:,n)= exp(real( repmat(z_G,1,Np).* repmat(exp(-1j*pp'),Ll,1)   ))' * (factG./(2*pi*besseli(0,abs(z_G)) ) );% + fact_inf./(2*pi); z=0 --> pdf=1./(2*pi)
    % %         beta_tkh(:,n)=beta_tkh(:,n)./sum(beta_tkh(:,n));
    % %
    % %
    % %
    % %         plot(pp,alpha_plot(:,n), 'r');
    % %         hold on;
    % %         plot(pp,alpha_exact(:,n), 'g');
    % %         plot(pp,alpha_tkh(:,n), 'k');
    % %         xlabel({'n= ',n});
    % %         hold off
    % %     legend('Tikh','exact','Gauss','Tikh-G')
    %%%_____________ END plot BETA __________
    
end
Papo=beta_plot.*alpha_plot;
Papo=Papo./repmat(sum(Papo,1),Np,1);

Papo_exact=beta_exact.*alpha_exact;
Papo_exact=Papo_exact./repmat(sum(Papo_exact,1),Np,1);

BCJR.Alpha_Tracking=alpha_plot;
BCJR.Beta_Tracking=beta_plot;


L_beta(pilot,:)=[];
L_alpha(pilot,:)=[];
Z_beta(pilot,:)=[];
Z_alpha(pilot,:)=[];
Fact_beta(pilot,:)=[];
Fact_alpha(pilot,:)=[];

% % Z_betaG(pilot,:)=[];
% % Z_alphaG(pilot,:)=[];
% % Fact_betaG(pilot,:)=[];
% % Fact_alphaG(pilot,:)=[];

% % BCJR.Z_betaG=Z_betaG;
% % BCJR.Fact_betaG=Fact_betaG;
% % BCJR.Z_alphaG=Z_alphaG;
% % BCJR.Fact_alphaG=Fact_alphaG;

BCJR.L_beta=L_beta;
BCJR.Z_beta=Z_beta;
BCJR.Fact_beta=Fact_beta;

BCJR.L_alpha=L_alpha;
BCJR.Z_alpha=Z_alpha;
BCJR.Fact_alpha=Fact_alpha;
BCJR.phase_shift=phase_shift.';


end
function z=gamma_fct(sig,z)

z=z./(1+abs(z).*sig);

end
