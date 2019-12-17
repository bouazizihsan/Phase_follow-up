function [BCJR]=Phase_track_BCJR_Gauss_Interp( y, N0, sig2, X, pilot_symbols,pilot,Lframe,payload,Phase_noise)
ss           =size(y);
D            =ss(2);
Ns           =ss(1);
pilot_pos    = zeros(Ns,2);
pilot_pos(pilot,:)= pilot_symbols;
M            =length(X);
% N0           =N0.*eye(D);
sigw2        =sig2.*eye(D);
Np           = 2^9;
Npilots      = length(pilot_symbols);
Phase_limit  = pi;%3*sqrt(sig2*Lframe/4);  % we suppose the phase varies in the interval (-Phase_limit, Phase_limit)

dpp         = 2.*Phase_limit./Np;
p           = linspace(-Phase_limit,Phase_limit,Np).';

[M1,M2]=meshgrid(p);
M11=M1(:);
M22=M2(:);


pdf_Gauss=@(t,m,s2,d) exp(-0.5*sum(((t-m )'*inv(s2).*(t-m ).'),2) )./((2*pi).^(d/2)*det(s2).^(0.5));

pdf_Gauss_log=@(t,m,s2,d) ((-0.5*(t-m )'*inv(s2)*(t-m )) - log( (2*pi).^(d/2)*det(s2).^(0.5) ));


%processing

phase_shift=interp1(pilot', Phase_noise(pilot,:), (1:Ns)) ;
y=y.*exp(-1j.*phase_shift);

Theta_zero=zeros(M,ss(2),Ns);
gamma=zeros(Np^D,Ns);
for n=1:Ns
    auxi=[];
    aux1=[];
    aux_SNL=[];
    
    if( norm(pilot_pos(n,:))~=0)
        X_Hard=pilot_pos(n,:);
        
    else
        
        X_Hard=X;
        
    end
    
    MM=size(X_Hard,1);
    Moy{1,n} =angle(repmat(y(n,:),MM,1)./X_Hard);
    
    for i=1:MM
%         aux1=[aux1;diag(N0./(2*abs(X_Hard(i,:)).^2) )];
%         aux_SNL=[aux_SNL;diag(N0./(2*abs(X_Hard(i,:)).*abs(y(n,:))))] ;
%         auxi=[auxi;(log(1./prod(( abs(X_Hard(i,:)).*(sqrt(abs(X_Hard(i,:).*y(n,:)))) ))) + pdf_Gauss_log( abs(y(n,:)./(X_Hard(i,:))).',1,aux1((i-1)*D+1:i*D,:),D) )];
    

        aux1=[aux1;diag(N0./(2*abs(X_Hard(i,:)).^2) )];
        aux_SNL=[aux_SNL;diag(N0./(2*abs(X_Hard(i,:)).^2))] ;
        auxi=[auxi;(log(1./(prod( abs(X_Hard(i,:)).^2 ))) + (pdf_Gauss_log( abs(y(n,:)./(X_Hard(i,:))).',1,aux1((i-1)*D+1:i*D,:),D))')];
 
    end
    Theta_zero(:,:,n)=[Moy{1,n};zeros(M-size(Moy{1,n},1),ss(2))];
    Xci_p0{1,n}=auxi;
    
    Snl{1,n}=aux_SNL;
    SigmaT{1,n}=aux1;
end


%ALPHA
alpha=zeros(Np^D,Ns);
alpha(:,1)=1;
alpha_p0{1,1}=[1 1] ;
alpha_p0{1,2}=diag(ones(1,D).*Inf);
% alpha_p0{1,3}=[1 1];

for n=2:Ns
    
   [Mnl,S]=Projection_of_productD(Moy{1,n-1},Snl{1,n-1},Xci_p0{1,n-1},alpha_p0{n-1,1},alpha_p0{n-1,2},sigw2,D,M1,M2,Phase_noise(n,:));
   
    alpha_p0{n,1}= Mnl;
    alpha_p0{n,2}=S; 
end

%beta
beta=zeros(Np^D,Ns);
beta(:,Ns)=1;
beta_p0{Ns,1}=[1 1] ;
beta_p0{Ns,2}=diag(ones(1,D).*Inf);

for n=Ns-1:-1:1
    [Mnl,S]=Projection_of_productD(Moy{1,n+1},Snl{1,n+1},Xci_p0{1,n+1},beta_p0{n+1,1},beta_p0{n+1,2},sigw2,D,M1,M2,Phase_noise(n,:));

    beta_p0{n,1}= Mnl;
    beta_p0{n,2}=S; 
      
end

%Product of Alpha and Beta ==> one gaussian

M_alpha_beta=zeros(Ns,D);
V_alpha_beta=zeros(Ns*D,D);

for i=1:Ns

    if(i==1)
        [Mnl,S]=Projection_of_productD(beta_p0{i,1},beta_p0{i,2},[0 0],alpha_p0{i,1},alpha_p0{i,2},0,D,M1,M2,Phase_noise(i,:));
        
    elseif(i==Ns)
        [Mnl,S]=Projection_of_productD(alpha_p0{i,1},alpha_p0{i,2},[0 0],beta_p0{i,1},beta_p0{i,2},0,D,M1,M2,Phase_noise(i,:));

    else % cas ou il y a pas de varience inf
        [Mnl,S]=Projection_of_productD(alpha_p0{i,1},alpha_p0{i,2},[0 0],beta_p0{i,1},beta_p0{i,2},0,D,M1,M2,Phase_noise(i,:));
    end
            
        V_alpha_beta((i-1)*D+1:i*D,:)= diag(diag(S));
        M_alpha_beta(i,:)= Mnl;
    
end

q=[payload*D-1 payload*D].';
q=q(:);
V_alpha_beta=V_alpha_beta(q,:);

BCJR.Mnl=Theta_zero(:,:,payload);
BCJR.Pmean=M_alpha_beta(payload,:);
BCJR.Pvar=V_alpha_beta;
BCJR.phase_shift=phase_shift;

BCJR.Snl=SigmaT;

BCJR.Pmean_apos1=M_alpha_beta(:,1);

BCJR.Pmean_apos2=M_alpha_beta(:,2);
end
