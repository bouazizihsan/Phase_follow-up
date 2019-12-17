function [BCJR]=DP( y, N0, sigw2, X, pilot_symbols,pilot)


Ns           =length(y);
pilot_pos    = zeros(Ns,1);
pilot_pos(pilot)= pilot_symbols;
M            =length(X);
Np           = 2^10;
%Npilots      = length(pilot_symbols);
Phase_limit  = pi;
dpp          = 2.*Phase_limit./Np;
p           = [linspace(-Phase_limit,0,Np/2) linspace(dpp,Phase_limit,Np/2)]';

pdf_Gauss   = @(t,m,s2) exp( -abs(t-m ).^2./(2*s2) )./sqrt(2*pi*s2); % pdf(t|m), normpdf(x,mu,sigma.^2)
g = fftshift(pdf_Gauss( p ,0 ,sigw2)); %repmat(-dphase_shift(n-1),Np,1 )


% yy=repmat(y.',Np,1);
% pp=repmat(p,1,Ns);
% YY=yy.*exp(-1j.*pp);
% 
% [XX,i]=Hard_Decision(YY(:),X); %% Hard_Decision(YY(:),X);  hard_QAM(YY(:),X);
% %[Z,J]=hard_QAM(YY(:),X);
% % ZZ=reshape(Z,Np,Ns);
% % JJ=reshape(J,Np,Ns);
% 
% XX=reshape(XX,Np,Ns);
% ii=reshape(i,Np,Ns);

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
        %         [cel{1,n},uniq]=unique(ii(:,n));
        X_Hard=X; % X(ii(uniq,n));
    end
    MM=length(X_Hard);
    for xh=1:MM
        gamma_exactX(:,n,xh) = pdf_Gauss( X_Hard(xh)*exp(1j*p) ,y(n), N0/2);
        gamma_exact(:,n) = gamma_exact(:,n) + gamma_exactX(:,n,xh);
    end
    gamma_exactX(:,n,:)=gamma_exactX(:,n,:)./sum(gamma_exact(:,n));
end


%ALPHA
for n=2:Ns
    
    %% calcul exact
    alpha_exact(:,n) = ifft( fft( alpha_exact(:,n-1).* gamma_exact(:,n-1), Np ) .* fft(g, Np)  );
    alpha_exact(:,n)=abs(alpha_exact(:,n)); %.*(alpha_exact(:,n)>0);
    alpha_exact(:,n) = alpha_exact(:,n)./sum( alpha_exact(:,n));
end

%beta
for n=Ns-1:-1:1
    
    
    %% calcul exact   
    beta_exact(:,n) = ifft( fft( beta_exact(:,n+1).* gamma_exact(:,n+1), Np ) .* fft(g, Np)  );
    beta_exact(:,n) = abs(beta_exact(:,n)); %.*(beta_exact(:,n)>0);
    %%normalize
    beta_exact(:,n) = beta_exact(:,n)./sum( beta_exact(:,n));
    
end

p_apo=beta_exact.*alpha_exact;
p_apo=p_apo./repmat(sum(p_apo,1),Np,1);

p_apoX=repmat(p_apo,1,1,M).*gamma_exactX;
p_apoX= p_apoX./repmat(sum(sum(p_apoX,3),1),Np,1,M);
p_apoXX(1:Ns,1:M)=sum(p_apoX,1);
p_apoXX=abs(p_apoXX); %.*(p_apoXX>0);

Pmean = ( p'*p_apo ); %*dpp;
Pvar  = ( (p.^2)'*p_apo )- Pmean.^2; %*dpp - Pmean.^2;
Pmean=Pmean.';
Pvar=Pvar.';

% Pmean(pilot)=[];
% Pvar(pilot)=[];
BCJR.p_apoX=p_apoXX;
clear -regexp  ^p_apoXX ^p_apoX ^beta_exact ^alpha_exact ^gamma_exact
BCJR.Pmean=Pmean;
BCJR.Pvar=Pvar;
end
