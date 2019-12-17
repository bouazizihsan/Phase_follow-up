function [BCJR]=Phase_track_Tikhonov( y, N0, sigw2, x, pilot_symbols,pilot,Lframe,NClass )

epsilon=2;
Ns           =length(y);
pilot_pos(pilot)= pilot_symbols;
M            =length(x);
Len            = NClass-1;
Np           = 2^9;
% BCJR.alpha=zeros(Np,Ns);
Npilots      = length(pilot_symbols);
Phase_limit  = 3*sqrt(sigw2*Lframe/4);  % we suppose the phase varies in the interval (-Phase_limit, Phase_limit)
dpp          = 2.*Phase_limit./Np;
pp           =   linspace(-Phase_limit,Phase_limit,Np)';
% ppp=[pp;pp+2.*pi];

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



alpha_plot=zeros(Np,Ns);
alpha_exact=zeros(Np,Ns);
alpha_exact(:,1)=1/(2*pi);
alpha_exact(:,1)=alpha_exact(:,1)./sum(alpha_exact(:,1));

beta_plot=zeros(Np,Ns);
beta_exact=zeros(Np,Ns);
beta_exact(:,Ns)=1/(2*pi);
beta_exact(:,Ns)=beta_exact(:,Ns)./sum(beta_exact(:,Ns));

gamma=zeros(Np,Ns);

for n=1:Ns
    if(pilot_pos(n)~=0)
        gamma(:,n) = pdf_Gauss( y(n) ,pilot_pos(n).*exp(1j*pp), N0/2);
    else
        for nn=1:M
            gamma(:,n) = gamma(:,n) + pdf_Gauss( y(n) , x(nn).*exp(1j*pp), N0/2);
        end
    end
    gamma(:,n)=gamma(:,n)./sum(gamma(:,n));
end

alpha=[1 0];
Z_alpha=zeros(Ns,NClass);
Fact_alpha=zeros(Ns,NClass);
Fact_alpha(1)=1;
% % 
% % alphaG=[1 0];
% % Z_alphaG=zeros(Ns,NClass);
% % Fact_alphaG=zeros(Ns,NClass);
% % Fact_alphaG(1)=1;

L_alpha=zeros(Ns,1);
L_alpha(1)=1;

L=size(alpha,1); % alpha  class size
Ll=1; % alpha exacte size

for n=2:Ns
    
    bessel2=besseli(0, abs(alpha(:,2)) );
    aux= find(bessel2 ==inf);
    if(aux)
        bessel2(aux)=exp(abs(alpha(aux,2)))./sqrt(2*pi*abs(alpha(aux,2)));
    end
    
    if(pilot_pos(n-1)~=0)
        z=(alpha(:,2)) + repmat((2/N0)*y(n-1)*conj(pilot_pos(n-1)),L,1);
        bessel1=besseli(0,abs(z));
        [l,c,v]= find(bessel1 ==inf);
        if(v)
            bessel1(l,c)=exp(abs(z(l,c)))./sqrt(2*pi*abs(z(l,c)));
        end
        
        lambda=  bessel1./bessel2.* exp( -abs(pilot_pos(n-1)).^2 /N0);
        lambda=(alpha(:,1)).*lambda;
    else
        z=repmat((alpha(:,2)),1,M) + repmat((2/N0)*y(n-1)*conj(x),L,1);
        bessel1=besseli(0,abs(z));
        [l,c,v]= find(bessel1 ==inf);
        if(v)
            bessel1(l,c)=exp(abs(z(l,c)))./sqrt(2*pi*abs(z(l,c)));
        end
        
        lambda= (1/M) .* bessel1./repmat(bessel2,1,M).* repmat(exp(-abs(x).^2/N0),size(z,1),1);
        lambda=repmat((alpha(:,1)),1,M).*lambda;
    end
    
    %___________________________
    
% %     
% %     bessel2=besseli(0, abs(alphaG(:,2)) );
% %     aux= find(bessel2 ==inf);
% %     if(aux)
% %         bessel2(aux)=exp(abs(alphaG(aux,2)))./sqrt(2*pi*abs(alphaG(aux,2)));
% %     end
% %     
% %     if(pilot_pos(n-1)~=0)
% %         zG=(alphaG(:,2)) + repmat((2/N0)*y(n-1)*conj(pilot_pos(n-1)),Ll,1);
% %         bessel1=besseli(0,abs(zG));
% %         [l,c,v]= find(bessel1 ==inf);
% %         if(v)
% %             bessel1(l,c)=exp(abs(zG(l,c)))./sqrt(2*pi*abs(zG(l,c)));
% %         end
% %         
% %         lambdaG=  bessel1./bessel2.* exp( -abs(pilot_pos(n-1)).^2 /N0);
% %         lambdaG=(alphaG(:,1)).*lambdaG;
% %     else
% %         zG=repmat((alphaG(:,2)),1,M) + repmat((2/N0)*y(n-1)*conj(x),Ll,1);
% %         bessel1=besseli(0,abs(zG));
% %         [l,c,v]= find(bessel1 ==inf);
% %         if(v)
% %             bessel1(l,c)=exp(abs(zG(l,c)))./sqrt(2*pi*abs(zG(l,c)));
% %         end
% %         
% %         lambdaG= (1/M) .* bessel1./repmat(bessel2,1,M).* repmat(exp(-abs(x).^2/N0),size(zG,1),1);
% %         lambdaG=repmat((alphaG(:,1)),1,M).*lambdaG;
% %     end
% %     zG =gamma_fct(sigw2,zG);
% %     bbb=find(lambdaG==inf);
    %___________________________
    z =gamma_fct(sigw2,z);
    [fact,zT]=KLCluster_U(lambda,z,epsilon,Len);
        % [fact,factG,zT,z_G]=KLCluster_U_new(lambda,lambdaG,z,zG,epsilon,Len);
    
    
    alpha = [fact zT];
%     alphaG=[factG z_G];
    L = size(alpha,1);
%     Ll=length(z_G);
    
    Z_alpha(n,1:L)=zT;
    Fact_alpha(n,1:L)=fact;
    
% %     Z_alphaG(n,1:Ll)=z_G;
% %     Fact_alphaG(n,1:Ll)=z_G;
    
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
L=size(beta,1);
Z_beta=zeros(Ns,NClass);
Fact_beta=zeros(Ns,NClass);
Fact_beta(Ns)=1;

betaG=[1 0];
Ll=size(beta,1);
Z_betaG=zeros(Ns,NClass);
Fact_betaG=zeros(Ns,NClass);
Fact_betaG(Ns)=1;

L_beta=zeros(Ns,1);
L_beta(Ns)=1;

for n=Ns-1:-1:1
    
    lambda=zeros(L,M);
    phi=zeros(L,M);
    z=zeros(L,M);
    
    bessel2=besseli(0, abs(beta(:,2)) );
    aux= find(bessel2 ==inf);
    if(aux)
        bessel2(aux)=exp(abs(beta(aux,2)))./sqrt(2*pi*abs(beta(aux,2)));
    end
    
    if(pilot_pos(n+1) ~= 0)
        z=(beta(:,2)) + repmat((2/N0)*y(n+1)*conj(pilot_pos(n+1)),L,1);
        bessel1=besseli(0,abs(z));
        [l,c,v]= find(bessel1 ==inf);
        if(v)
            bessel1(l,c)=exp(abs(z(l,v)))./sqrt(2*pi*abs(z(l,v)));
        end
        %bessel1=reshape(bessel1,size(z,1),[]);
        
        lambda=  bessel1./bessel2.* exp(-abs(pilot_pos(n+1)).^2/N0);
        lambda=(beta(:,1)).*lambda;
        
    else
        
        z=repmat((beta(:,2)),1,M) + repmat((2/N0)*y(n+1)*conj(x),L,1);
        bessel1=besseli(0,abs(z));
        [l,c,v]= find(bessel1 ==inf);
        if(v)
            bessel1(l,c)=exp(abs(z(l,v)))./sqrt(2*pi*abs(z(l,v)));
        end
        
        lambda= (1/M) .* bessel1./repmat(bessel2,1,M) .* repmat(exp(-abs(x).^2/N0),size(z,1),1);
        lambda=repmat((beta(:,1)),1,M).*lambda;
    end
    
    %___________________________
    
% %     
% %     bessel2=besseli(0, abs(betaG(:,2)) );
% %     aux= find(bessel2 ==inf);
% %     if(aux)
% %         bessel2(aux)=exp(abs(betaG(aux,2)))./sqrt(2*pi*abs(betaG(aux,2)));
% %     end
% %     
% %     if(pilot_pos(n+1)~=0)
% %         zG=(betaG(:,2)) + repmat( (2/N0)*y(n+1)*conj(pilot_pos(n+1)),Ll,1);
% %         bessel1=besseli(0,abs(zG));
% %         [l,c,v]= find(bessel1 ==inf);
% %         if(v)
% %             bessel1(l,c)=exp(abs(zG(l,c)))./sqrt(2*pi*abs(zG(l,c)));
% %         end
% %         
% %         lambdaG=  bessel1./bessel2.* exp( -abs(pilot_pos(n+1)).^2 /N0);
% %         if( size(alphaG(:,1))~= size(lambdaG) )
% %             stop=1;
% %         end
% %         lambdaG=(betaG(:,1)).*lambdaG;
% %     else
% %         
% %         zG=repmat((betaG(:,2)),1,M) + repmat((2/N0)*y(n+1)*conj(x),Ll,1);
% %         bessel1=besseli(0,abs(zG));
% %         [l,c,v]= find(bessel1 ==inf);
% %         if(v)
% %             bessel1(l,c)=exp(abs(zG(l,c)))./sqrt(2*pi*abs(zG(l,c)));
% %         end
% %         
% %         lambdaG= (1/M) .* bessel1./repmat(bessel2,1,M).* repmat(exp(-abs(x).^2/N0),size(zG,1),1);
% %         lambdaG=repmat((betaG(:,1)),1,M).*lambdaG;
% %     end
% %     zG=gamma_fct(sigw2,zG);
% %     bbb=find(lambdaG==inf);
    %___________________________
    z=gamma_fct(sigw2,z);
    
    aaa=find(lambda==inf);
    
    if(aaa)
        stop=1;
% %     elseif(bbb)
% %         stop=1;
    
    else
        
        [fact,zT]=KLCluster_U(lambda,z,epsilon,Len);
% %         [fact,factG,zT,z_G]=KLCluster_U_new(lambda,lambdaG,z,zG,epsilon,Len);
    end
    beta=[fact zT];
% %     betaG=[factG z_G];    
% %     Ll=size(betaG,1);
    
    L=size(beta,1);
    Z_beta(n,1:L)=zT;
    Fact_beta(n,1:L)=fact;
    L_beta(n)=L;
    
% %     Z_betaG(n,1:Ll)=z_G;
% %     Fact_betaG(n,1:Ll)=factG;
    
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
