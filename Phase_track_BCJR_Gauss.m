function [BCJR]=Phase_track_BCJR_Gauss( y, N0, sigw2, X, pilot_symbols,pilot,Lframe,Gamma_Apprx,Pilot_Gamma_Apprx,AlphaGammaProj,Phase_noise)
Ns           =length(y);
pilot_pos    = zeros(Ns,1);
pilot_pos(pilot)= pilot_symbols;
M            =length(X);

NP=length(pilot);
[BCJRP]=BCJR_Pilot( y(pilot), N0, Lframe*sigw2, pilot_symbols,Pilot_Gamma_Apprx);

Np           = 2^9;
%Npilots      = length(pilot_symbols);
Phase_limit  = pi;
dpp          = 2.*Phase_limit./Np;
p           = [linspace(-Phase_limit,0,Np/2) linspace(dpp,Phase_limit,Np/2)]';

pdf_Gauss   = @(t,m,s2) exp( -abs(t-m ).^2./(2*s2) )./sqrt(2*pi*s2); % pdf(t|m), normpdf(x,mu,sigma.^2)
g = fftshift(pdf_Gauss( p ,0 ,sigw2)); %repmat(-dphase_shift(n-1),Np,1 )

%processing


new_phases=BCJRP.apo(:,1); % unwrap(angle(y(pilot,:)./pilot_symbols)); %   Phase_noise(pilot); %
phase_shift=(interp1(pilot', new_phases, (1:Ns))).' ;
% phase_shift=phase_shift*0;
% phase_shift=Phase_noise;
dd_phase_shift = diff( phase_shift );
% y=y.*exp(-1j*phase_shift);
phase_ref=(phase_shift);

gamma_exact = zeros(Np,Ns);
gamma_exactX = zeros(Np,Ns,M);
alpha_exact=zeros(Np,Ns);
beta_exact=zeros(Np,Ns);

alpha_exact(:,1)=1;
beta_exact(:,Ns)=1;


for n=1:Ns
    
    if( pilot_pos(n)~=0)
        X_Hard=pilot_pos(n);
        
    else
        
        X_Hard=X;
        
    end
    MM=length(X_Hard);
    for xh=1:MM
        gamma_exactX(:,n,xh) = pdf_Gauss( X_Hard(xh)*exp(1j*p) ,y(n), N0/2);
        gamma_exact(:,n) = gamma_exact(:,n) + gamma_exactX(:,n,xh);
    end
    gamma_exactX(:,n,:)=gamma_exactX(:,n,:)./sum(gamma_exact(:,n));
    
    
    switch Gamma_Apprx
        case 'SP'
%             m=unwrap([repmat(phase_shift(n),1,length(X_Hard));(angle(y(n)./(X_Hard))).']);
%             Mnl_p0{1,n} =m(2,:).';
            Mnl_p0{1,n} =angle(y(n)./(X_Hard));
            Snl{1,n} = N0./(2*abs(X_Hard).*abs(y(n)));
            Xci_p0{1,n}=pdf_Gauss( abs(y(n)) ,abs(X_Hard),N0/2)./ (sqrt(abs(X_Hard*y(n))) );
        case 'LTmax'
%             m=unwrap([repmat(phase_shift(n),1,length(X_Hard));(angle(y(n)./(X_Hard))).']);
%             Mnl_p0{1,n} =m(2,:).';
            Mnl_p0{1,n} =angle(y(n)./(X_Hard));
            Snl{1,n} = N0./(2*abs(X_Hard).^2);
            Xci_p0{1,n}=pdf_Gauss( abs(y(n)) ,abs(X_Hard),N0/2)./ (sqrt(abs(X_Hard*y(n))) );
        case 'LT'
            z=conj(y(n))*exp(1j*phase_ref(n));
            Mnl_p0{1,n} =phase_ref(n)-imag(z.*X_Hard)./(abs(X_Hard).^2);
            Snl{1,n} =  N0./(2.*abs(X_Hard).^2);
            Xci_p0{1,n}=pdf_Gauss( real( z*X_Hard ) ,abs(X_Hard).^2,  N0./2.*(abs(X_Hard).^2));
            Xci_p0{1,n}=Xci_p0{1,n}./sum(Xci_p0{1,n});
        otherwise
            error('unknown data Gamma Apprx method ');
            
    end
    
    
    %     if(~(is_BCJRP))
    %         Gamma_Plot(X_Hard,(Xci_p0{1,n}./sum(Xci_p0{1,n})),Mnl_p0{1,n},Snl{1,n},'b')
    %     end
    
    
end
phase_shift=phase_shift*0;
dd_phase_shift = diff( phase_shift );
phase_ref=(phase_shift-phase_shift);

%ALPHA
parameters1{1,1}=0;
parameters1{1,2}=inf;
parameters1{1,3}=log(1);

parameters1{Ns,1}=0;
parameters1{Ns,2}=0;
parameters1{Ns,3}=0;
for n=2:Ns-1
    switch AlphaGammaProj
        case 'OneSideCorr'
      
            nn=min(find(pilot>n));
%             m=unwrap([phase_shift(n) BCJRP.beta_gamma(nn,1)]);
%             Corr=[-diff(m);BCJRP.beta_gamma(nn,2)+(pilot(nn)-n).*sigw2;0]; % (nn-1)*Lframe+1 = pilot(nn), (pilot>=n) correction not affect pilot position, (pilot>n) correction affect pilot position
            Corr=[BCJRP.beta_gamma(nn,1)-phase_shift(n);BCJRP.beta_gamma(nn,2)+(pilot(nn)-n).*sigw2;0]; % (nn-1)*Lframe+1 = pilot(nn), (pilot>=n) correction not affect pilot position, (pilot>n) correction affect pilot position
                         
        case 'BothSideCorr'
            
            nn=min(find(pilot>n));
%             m=unwrap([phase_shift(n) BCJRP.beta_gamma(nn,1)]);
%             beta_pilot=[-diff(m);BCJRP.beta_gamma(nn,2)+(pilot(nn)-n).*sigw2;0];
            beta_pilot=[BCJRP.beta_gamma(nn,1)-phase_shift(n);BCJRP.beta_gamma(nn,2)+(pilot(nn)-n).*sigw2;0];
%             m=unwrap([phase_shift(n) BCJRP.alpha_gamma(nn-1,1)]);
%             alpha_pilot=[-diff(m);BCJRP.alpha_gamma(nn-1,2)+(n-pilot(nn-1)).*sigw2;0];
            alpha_pilot=[BCJRP.alpha_gamma(nn-1,1)-phase_shift(n);BCJRP.alpha_gamma(nn-1,2)+(n-pilot(nn-1)).*sigw2;0];
            [Corr(1),Corr(2),Corr(3)]=Produit_Alpha_Gamma(alpha_pilot(1) , alpha_pilot(2) ,alpha_pilot(3),beta_pilot(1),beta_pilot(2),beta_pilot(3));
        otherwise
            Corr=0;
    end
    
    [Mnl,Sigma,ph]=ProductProjection(Mnl_p0{1,n-1},Snl{1,n-1},log(Xci_p0{1,n-1}),parameters1{n-1,1},parameters1{n-1,2},parameters1{n-1,3},Corr,AlphaGammaProj);
    
    parameters1{n,1}=Mnl-dd_phase_shift(n-1);
    parameters1{n,2}=Sigma+sigw2;
    parameters1{n,3}=ph;
    
%     Alpha=pdf_Gauss(p,parameters1{n,1},parameters1{n,2});
%     Alpha=Alpha./sum(Alpha);

    vars = {'nn','m','Corr','ph','Mnl','Sigma','ph'};
    clear(vars{:})
    
    %% calcul exact
    alpha_exact(:,n) = ifft( fft( alpha_exact(:,n-1).* gamma_exact(:,n-1), Np ) .* fft(g, Np)  );
    alpha_exact(:,n)=abs(alpha_exact(:,n)); %.*(alpha_exact(:,n)>0);
    alpha_exact(:,n) = alpha_exact(:,n)./sum( alpha_exact(:,n));
    
%     figure(1)
%     plot(p,alpha_exact(:,n),'g')
%     hold on
%     plot(p,Alpha,'r')
%     hold off
end

%beta
parameters2{Ns,1}=0;
parameters2{Ns,2}=inf;
parameters2{Ns,3}=log(1);

parameters2{1,1}=0;
parameters2{1,2}=0;
parameters2{1,3}=0;
for n=Ns-1:-1:2
    switch AlphaGammaProj
        case 'OneSideCorr'           
            nn=max(find(pilot<n));
%             m=unwrap([phase_shift(n) BCJRP.alpha_gamma(nn,1)]);
%             Corr=[-diff(m);BCJRP.alpha_gamma(nn,2)+(n-pilot(nn)).*sigw2;0];
            
            Corr=[BCJRP.alpha_gamma(nn,1)-phase_shift(n);BCJRP.alpha_gamma(nn,2)+(n-pilot(nn)).*sigw2;0];
        case 'BothSideCorr'
            
            nn=max(find(pilot<n));
%             m=unwrap([phase_shift(n) BCJRP.alpha_gamma(nn,1)]);
%             alpha_pilot=[-diff(m);BCJRP.alpha_gamma(nn,2)+(n-pilot(nn)).*sigw2;0];
            alpha_pilot=[BCJRP.alpha_gamma(nn,1)-phase_shift(n);BCJRP.alpha_gamma(nn,2)+(n-pilot(nn)).*sigw2;0];
%             m=unwrap([phase_shift(n) BCJRP.beta_gamma(nn+1,1)]);
%             beta_pilot=[-diff(m);BCJRP.beta_gamma(nn+1,2)+(pilot(nn+1)-n).*sigw2;0];
            beta_pilot=[BCJRP.beta_gamma(nn+1,1)-phase_shift(n);BCJRP.beta_gamma(nn+1,2)+(pilot(nn+1)-n).*sigw2;0];
            [Corr(1),Corr(2),Corr(3)]=Produit_Alpha_Gamma(alpha_pilot(1) , alpha_pilot(2) ,alpha_pilot(3),beta_pilot(1),beta_pilot(2),beta_pilot(3));
        otherwise
            Corr=0;
    end
    
    [Mnl,Sigma,ph]=ProductProjection(Mnl_p0{1,n+1},Snl{1,n+1},log(Xci_p0{1,n+1}),parameters2{n+1,1},parameters2{n+1,2},parameters2{n+1,3},Corr,AlphaGammaProj);
    parameters2{n,1}=Mnl+dd_phase_shift(n);
    parameters2{n,2}=Sigma+sigw2;
    parameters2{n,3}=ph;
    
%     Beta=pdf_Gauss(p,parameters1{n,1},parameters1{n,2});
%     Beta=Beta./sum(Beta);
    
    vars = {'nn','m','Corr','ph','Mnl','Sigma','ph'};
    clear(vars{:})
    
     %% calcul exact   
    beta_exact(:,n) = ifft( fft( beta_exact(:,n+1).* gamma_exact(:,n+1), Np ) .* fft(g, Np)  );
    beta_exact(:,n) = abs(beta_exact(:,n)); %.*(beta_exact(:,n)>0);
    %%normalize
    beta_exact(:,n) = beta_exact(:,n)./sum( beta_exact(:,n));
    
%     figure(1)
%     plot(p,beta_exact(:,n),'g')
%     hold on
%     plot(p,Beta,'r')
%     hold off
    
end

Pmean=zeros(Ns,1);
Pvar=zeros(Ns,1);

for n=1:Ns
    if(n==Ns)
        [Pmean(n,1),Pvar(n,1),ph]=ProductProjection(parameters1{n,1},parameters1{n,2},parameters1{n,3},parameters2{n,1},parameters2{n,2},parameters2{n,3},0,'None');
    else
        [Pmean(n,1),Pvar(n,1),ph]=ProductProjection(parameters2{n,1},parameters2{n,2},parameters2{n,3},parameters1{n,1},parameters1{n,2},parameters1{n,3},0,'None');
    end
    
end

p_apo=beta_exact.*alpha_exact;
p_apo=p_apo./repmat(sum(p_apo,1),Np,1); % /dpp;
p_apoX=repmat(p_apo,1,1,M).*gamma_exactX;
p_apoX= p_apoX./repmat(sum(sum(p_apoX,3),1),Np,1,M);
for n=1:Ns
AlphaBeta=pdf_Gauss(p,Pmean(n),Pvar(n,1));  
AlphaBeta=AlphaBeta./sum(AlphaBeta);
var=sum(p_apoX(:,n,:),3);
plot(p,var,'g')
hold on
plot(p,AlphaBeta,'r')
hold off
end
p_apoXX(1:Ns,1:M)=sum(p_apoX,1);
p_apoXX=abs(p_apoXX); %.*(p_apoXX>0);


BCJR.Pmean=Pmean;
BCJR.Pvar=Pvar;
BCJR.phase_shift=phase_shift;
clear -regexp  ^Pmean ^Pvar ^paramet
end
