function [BCJR]=Phase_track_BCJR_Gauss_App1( y, N0, sigw2, X, pilot_symbols,pilot,Lframe,Phase_noise ,hMap,Data,hDec,Nc,Nb)

% vth      = 6; % taux d'erreur= 25% eps=0.25; pour eps=0.3(30%) on a vth = 9.2; v=(-pi./(qfuncinv(1-eps/2)))^2
Ns           =length(y);
pilot_pos    = zeros(Ns,1);
pilot_pos(pilot)= pilot_symbols;
M            =length(X);
Np           = 2^9;
%Npilots      = length(pilot_symbols);
Phase_limit  =pi;%3*sqrt(sigw2*Lframe/4);  % we suppose the phase varies in the interval (-Phase_limit, Phase_limit)

dpp          = 2.*Phase_limit./Np;
p           = linspace(-Phase_limit,Phase_limit,Np)';
% p=[p-2*pi;p;p+2*pi];
% Np=Np*3;
pdf_Gauss   = @(t,m,s2) exp( -abs(t-m ).^2./(2*s2) )./sqrt(2*pi*s2); % pdf(t|m), normpdf(x,mu,sigma.^2)
pdf_Gauss_log   = @(t,m,s2) ( -abs(t-m ).^2./(2*s2) )-0.5*log((2*pi*s2));
pfd_xy=@(t,m,p,s2) exp( -abs(t-m.*exp(1j*p)).^2./(pi*s2) );
g = fftshift(pdf_Gauss( p ,0 ,sigw2));
gamma_exact = zeros(Np,Ns);
gamma_exactX= zeros(Np,Ns,M);
alpha_exact=zeros(Np,Ns);
beta_exact=zeros(Np,Ns);
alpha_exact(:,1)=1;
beta_exact(:,Ns)=1;

CC={'r--','b--','g--','k--'};
CCC={'r','b','k','g'};

%processing
phase=angle(y(pilot,:)./pilot_symbols);
new_phases=unwrap(phase);
phase_shift=(interp1(pilot', new_phases, (1:Ns))).' ;
dd_phase_shift = diff( phase_shift );
Y=y;
y=y.*exp(-1j*phase_shift);
% XX=repmat(X.',Ns,1);
% YY=repmat(y,1,M);
% Moy= angle(YY./(XX));
% Moy=unwrap(Moy);
% Moy(pilot,:)=repmat(new_phases,1,M);
phase_ref=(phase_shift-phase_shift);%phase_shift %
for n=1:Ns
    z=conj(y(n))*exp(1j*phase_ref(n));
    if( pilot_pos(n)~=0)
        X_Hard=pilot_pos(n);
%         Mnl_p0{1,n} =0; 
    else
        
        X_Hard=X;
%         Mnl_p0{1,n} = angle(y(n)./(X_Hard)); %imag(y(n)./(X_Hard));%
    end
    
     Mnl_p0{1,n} = angle(Y(n)./(X_Hard));%(178*pi/180); %imag(y(n)./(X_Hard));%
    %     %     Mnl_p0{1,n} =phase_shift(n)-imag(y(n).*X_Hard)./(abs(X_Hard).^2);
    %     %     Snl_Xci{1,n} = N0/2.*abs(X_Hard).^2;
    
    Snl_Xci{1,n} = N0./(2.*abs(X_Hard).^2);
    %     if((abs(X_Hard)>abs(y(n))) .* (pilot_pos(n)~=0))
    %         Snl{1,n} = N0./(2*abs(X_Hard).^2);
    %     else
    Snl{1,n} = N0./(2*abs(X_Hard).*abs(Y(n)));% Snl_Xci{1,n};%
    
    %     end

    Xci_p0{1,n}=pdf_Gauss( abs( Y(n)./X_Hard ) ,1, Snl_Xci{1,n})./ (sqrt(abs(X_Hard*Y(n))).*abs(X_Hard) ); % pdf_Gauss( real( y(n)./X_Hard ) ,1, Snl_Xci{1,n})./ (sqrt(abs(X_Hard*y(n))).*abs(X_Hard) );%
    TawX{1,n}=exp(-(1./N0)*(abs(Y(n)) + abs(X_Hard)).^2);
    Xci_p0{1,n}=Xci_p0{1,n}-TawX{1,n}./sqrt(abs(X_Hard).*abs(Y(n)));
    %     Mnl_p0_LT{1,n} =imag(y(n)./(X_Hard));%
    %     Snl_Xci_LT{1,n} =N0./(2.*abs(X_Hard).^2);
    %     Snl_LT{1,n} = Snl_Xci_LT{1,n};
    %     Xci_p0_LT{1,n}= pdf_Gauss( real( y(n)./X_Hard ) ,1, Snl_Xci_LT{1,n})./ (abs(X_Hard).^2 );%
    
    Mnl_p0_LT{1,n} =phase_ref(n)-imag(z.*(X_Hard))./abs(X_Hard).^2;
    Snl_Xci_LT{1,n} = N0/2.*abs(X_Hard).^2;
    Snl_LT{1,n} = N0./(2*abs(X_Hard).^2);
    Xci_p0_LT {1,n}= pdf_Gauss( real( z*X_Hard ) ,abs(X_Hard).^2, Snl_Xci{1,n});
    
    
    
    for nn=1:length(X_Hard)
        gamma_exactX(:,n,nn)= pdf_Gauss( X_Hard(nn)*exp(1j*p) ,Y(n), N0/2);
        gamma_exact(:,n) = gamma_exact(:,n) + gamma_exactX(:,n,nn);
        
    end
    gamma_exact(:,n)=gamma_exact(:,n)./sum(gamma_exact(:,n));
    
    gamma=(repmat(Xci_p0{1,n}.',Np,1)+pdf_Gauss_log(repmat(p,1,length(X_Hard)),repmat(Mnl_p0{1,n}.',Np,1),repmat(Snl{1,n}.',Np,1)));
    
    ggg=gamma(:);
    gamma=exp(gamma-jac_logH(ggg.'));
    Gamma{1,n}=gamma;
    
    gamma_LT=(repmat(Xci_p0_LT{1,n}.',Np,1)+pdf_Gauss_log(repmat(p,1,length(X_Hard)),repmat(Mnl_p0_LT{1,n}.',Np,1),repmat(Snl_LT{1,n}.',Np,1)));
    ggg=gamma_LT(:);
    gamma_LT=exp(gamma_LT-jac_logH(ggg.'));
    Sgamma_LT=sum(gamma_LT,2);
    Sgamma=sum(gamma,2);
    Sgamma_LT=Sgamma_LT./sum(Sgamma_LT);
    Sgamma=Sgamma./sum(Sgamma);
    
    f_gamma(n)=sum(Sgamma.*log(Sgamma./gamma_exact(:,n)));  %DKL(A,B,C,D);
    f_gamma_LT(n)=sum(Sgamma_LT.*log(Sgamma_LT./gamma_exact(:,n)));  %DKL(A,B,C,D);
    
    %     for i=1:length(X_Hard)
    %         hold on
    %         plot(p,gamma(:,i),CC{i});
    %         plot(p,gamma_LT(:,i),CCC{i});
    %     end
    %     Phase_noise(n)
    %     plot(Phase_noise(n),0,'r*');
    %     plot(p,gamma_exact(:,n),'m');
    %     title(n);
    %     clf
end

%ALPHA
parameters1{1,1}=0;
parameters1{1,2}=inf;
parameters1{1,3}=log(1);
parameters1{1,4}=0; % TawAlpha

parameters1_LT{1,1}=0;
parameters1_LT{1,2}=inf;
parameters1_LT{1,3}=log(1);

for n=2:Ns
%     if(pilot_pos(n-1)~=0)
%         Mnl=Mnl_p0{1,n-1};
%         Sigma=Snl{1,n-1};
%         ph=0;
%     else
    [Mnl,Sigma,ph,Taw]=Projection_of_product(Mnl_p0{1,n-1},Snl{1,n-1},log(Xci_p0{1,n-1}),TawX{1,n-1},parameters1{n-1,1},parameters1{n-1,2},parameters1{n-1,3},parameters1{n-1,4},2,sigw2,Phase_noise(n));
%     end
    parameters1{n,1}=Mnl; %-dd_phase_shift(n-1);
    parameters1{n,2}=Sigma;
    parameters1{n,3}=ph;
    parameters1{n,4}=Taw;
    [Mnl_LT,Sigma_LT,ph_LT,Taw]=Projection_of_product(Mnl_p0_LT{1,n-1},Snl_LT{1,n-1},log(Xci_p0_LT{1,n-1}),0,parameters1_LT{n-1,1},parameters1_LT{n-1,2},parameters1_LT{n-1,3},0,2,sigw2,Phase_noise(n));
    parameters1_LT{n,1}=Mnl_LT-dd_phase_shift(n-1);
    parameters1_LT{n,2}=Sigma_LT;
    parameters1_LT{n,3}=ph_LT;

    
    alpha_exact(:,n) = ifft( fft( alpha_exact(:,n-1).* gamma_exact(:,n-1), Np ) .* fft(g, Np)  );
    alpha_exact(:,n)=alpha_exact(:,n)./sum(alpha_exact(:,n));
    
    u_exact_alpha(:,n)= (cos(p)')*alpha_exact(:,n);
    var_exact_alpha(:,n)=((cos(p)').^2)*alpha_exact(:,n)-u_exact_alpha(:,n).^2;
%     u_exact_alpha(:,n)=(1-var_exact_alpha(:,n)/2).*exp(1j*u_exact_alpha(:,n));
    R=Operator(repmat(p,1,length(parameters1{n,3}))-repmat(parameters1{n,1}.',Np,1));
    alpha=(repmat(parameters1{n,3}.',Np,1)+pdf_Gauss_log(R,0,repmat(parameters1{n,2}.',Np,1)));
    ggg=alpha(:);
    alpha=exp(alpha-jac_logH(ggg.'));
    R=Operator(repmat(p,1,length(parameters1_LT{n,3}))-repmat(parameters1_LT{n,1}.',Np,1));
    alpha_LT=(repmat(parameters1_LT{n,3}.',Np,1)+pdf_Gauss_log(R,0,repmat(parameters1_LT{n,2}.',Np,1)));
    ggg=alpha_LT(:);
    alpha_LT=exp(alpha_LT-jac_logH(ggg.'));
    
    alpha_exact(:,n) =alpha_exact(:,n) .*(alpha_exact(:,n) >0);
    fi=find((alpha_exact(:,n)>0).*(alpha>(1e-308)));
    fi_LT=find((alpha_exact(:,n)>0).*(alpha_LT>(1e-308)));
    f_alpha(n)= sum(alpha_exact(fi,n).*log(alpha_exact(fi,n)./alpha(fi)));
    f_alpha_LT(n)=sum(alpha_exact(fi_LT,n).*log(alpha_exact(fi_LT,n)./alpha_LT(fi_LT)));
    if(n>400)
        figure(10)
        for i=1:length(parameters1{n,3})
            hold on
            plot(p,alpha(:,i),CC{i});
            plot(p,alpha_LT(:,i),CCC{i});
        end
        plot(p,alpha_exact(:,n),'g');
        plot(Phase_noise(n),0,'r*');
        title(n);
        
    end
end

%beta
parameters2{Ns,1}=0;
parameters2{Ns,2}=inf;
parameters2{Ns,3}=log(1);
parameters2{Ns,4}=0;

parameters2_LT{Ns,1}=0;
parameters2_LT{Ns,2}=inf;
parameters2_LT{Ns,3}=log(1);
u_exact_beta=zeros(1,Ns);
for n=Ns-1:-1:1
%     if(pilot_pos(n+1)~=0)
%         Mnl=Mnl_p0{1,n+1};
%         Sigma=Snl{1,n+1};
%         ph=0;
%     else
        [Mnl,Sigma,ph,Taw]=Projection_of_product(Mnl_p0{1,n+1},Snl{1,n+1},log(Xci_p0{1,n+1}),TawX{1,n+1},parameters2{n+1,1},parameters2{n+1,2},parameters2{n+1,3},parameters2{n+1,4},2,sigw2,Phase_noise(n));    
%     end
    parameters2{n,1}=Mnl; %+dd_phase_shift(n);
    parameters2{n,2}=Sigma;
    parameters2{n,3}=ph;
    parameters2{n,4}=Taw;
    [Mnl_LT,Sigma_LT,ph_LT,Taw]=Projection_of_product(Mnl_p0_LT{1,n+1},Snl_LT{1,n+1},log(Xci_p0_LT{1,n+1}),0,parameters2_LT{n+1,1},parameters2_LT{n+1,2},parameters2_LT{n+1,3},0,2,sigw2,Phase_noise(n));
    parameters2_LT{n,1}=Mnl_LT+dd_phase_shift(n);
    parameters2_LT{n,2}=Sigma_LT;
    parameters2_LT{n,3}=ph_LT;
    
    beta_exact(:,n) = ifft( fft( beta_exact(:,n+1).* gamma_exact(:,n+1), Np ) .* fft(g, Np)  );
    beta_exact(:,n)=beta_exact(:,n)./sum(beta_exact(:,n));
    
    u_exact_beta(:,n)= (cos(p)')*beta_exact(:,n);
    var_exact_beta(:,n)=((cos(p)').^2)*beta_exact(:,n)-u_exact_beta(:,n).^2;
%     u_exact_beta(:,n)=(1-var_exact_beta(:,n)/2).*exp(1j*u_exact_beta(:,n));
    R=Operator(repmat(p,1,length(parameters2{n,3}))-repmat(parameters2{n,1}.',Np,1));
    beta=(repmat(parameters2{n,3}.',Np,1)+pdf_Gauss_log(R,0,repmat(parameters2{n,2}.',Np,1)));
    ggg=beta(:);
    beta=exp(beta-jac_logH(ggg.'));
    R=Operator(repmat(p,1,length(parameters2_LT{n,3}))-repmat(parameters2_LT{n,1}.',Np,1));
    beta_LT=(repmat(parameters2_LT{n,3}.',Np,1)+pdf_Gauss_log(R,0,repmat(parameters2_LT{n,2}.',Np,1)));
    ggg=beta_LT(:);
    beta_LT=exp(beta_LT-jac_logH(ggg.'));
    beta_exact(:,n) =beta_exact(:,n) .*(beta_exact(:,n) >0);
    fi=find((beta_exact(:,n)>0).*(beta>(1e-308)));
    fi_LT=find((beta_exact(:,n)>0).*(beta_LT>(1e-308)));
    f_beta(n)= sum(beta_exact(fi,n).*log(beta_exact(fi,n)./beta(fi)));
    f_beta_LT(n)=sum(beta_exact(fi_LT,n).*log(beta_exact(fi_LT,n)./beta_LT(fi_LT)));
%     if(n<400)
%         figure(10)
%         for i=1:length(parameters2{n,3})
%             hold on
%             plot(p,beta(:,i),CC{i});
%             plot(p,beta_LT(:,i),CCC{i});
%     
%         end
%         plot(p,beta_exact(:,n),'g');
%         plot(Phase_noise(n),0,'r*');
%         title(n);
%         clf
%     end
%         figure(10)
%         for i=1:length(parameters2{n,3})
%             hold on
%             plot(p,beta(:,i),CC{i});
%             plot(p,beta_LT(:,i),CCC{i});
%     
%         end
%         plot(p,beta_exact(:,n),'g');
%         plot(Phase_noise(n),0,'r*');
%         title(n);
%         clf
end
alpha_beta_exact=beta_exact.*alpha_exact;
u_exact= (cos(p)')*alpha_beta_exact;
var_exact=((cos(p)').^2)*alpha_beta_exact-u_exact.^2;
% u_exact=(1-var_exact/2).*exp(1j*u_exact);

alpha_beta_exact=alpha_beta_exact.*(alpha_beta_exact>0);
alpha_beta_exact=alpha_beta_exact./repmat(sum(alpha_beta_exact),Np,1);
alpha_gamma_beta_exact(1:Ns,1:M)=sum(repmat(beta_exact,1,1,M).*gamma_exactX.*repmat(alpha_exact,1,1,M));

PmeanT=zeros(Ns,1);
PvarT=zeros(Ns,1);
Taw=zeros(Ns,1);

f_alpha_beta=zeros(1,Ns);
f_alpha_beta_LT=zeros(1,Ns);
for n=1:Ns
    
    if(n==Ns)
        [PmeanT(n,1),PvarT(n,1),ph,Taw(n,1)]=Projection_of_product(parameters1{n,1},parameters1{n,2},parameters1{n,3},parameters1{n,4},parameters2{n,1},parameters2{n,2},parameters2{n,3},parameters2{n,4},2,0,Phase_noise);
    else
        [PmeanT(n,1),PvarT(n,1),ph,Taw(n,1)]=Projection_of_product(parameters2{n,1},parameters2{n,2},parameters2{n,3},parameters2{n,4},parameters1{n,1},parameters1{n,2},parameters1{n,3},parameters1{n,4},2,0,Phase_noise);
    end
    
    if(n==Ns)
        [Pmean_LTT(n,:),Pvar_LTT(n,:),ph_LT,taux]=Projection_of_product(parameters1_LT{n,1},parameters1_LT{n,2},parameters1_LT{n,3},0,parameters2_LT{n,1},parameters2_LT{n,2},parameters2{n,3},0,2,0,Phase_noise);
    else
        [Pmean_LTT(n,:),Pvar_LTT(n,:),ph_LT,taux]=Projection_of_product(parameters2_LT{n,1},parameters2_LT{n,2},parameters2_LT{n,3},0,parameters1_LT{n,1},parameters1_LT{n,2},parameters1{n,3},0,2,0,Phase_noise);
    end
    R=Operator(p-PmeanT(n,:));
    alpha_beta(:,n)=pdf_Gauss(R,0,PvarT(n,:));
     R=Operator(p-Pmean_LTT(n,:));
    alpha_beta_LT(:,n)=pdf_Gauss(R,0,Pvar_LTT(n,:));
    
    alpha_beta(:,n)=alpha_beta(:,n)./sum(alpha_beta(:,n));
    alpha_beta_LT(:,n)=alpha_beta_LT(:,n)./sum(alpha_beta_LT(:,n));
    %     figure(2)
    %     plot(p,alpha_beta(:,n),CC{1});
    %
    %     hold on
    %     plot(p,alpha_beta_LT(:,n),CCC{1});
    %
    %     plot(p,alpha_beta_exact(:,n),'g');
    %     plot(Phase_noise(n),0,'r*');
    fi=find((alpha_beta_exact(:,n)>0).*(alpha_beta(:,n)>(1e-308)));
    fi_LT=find( (alpha_beta_exact(:,n)>0).*(alpha_beta_LT(:,n)>(1e-308)) );
    f_alpha_beta(n)=  sum(alpha_beta_exact(fi,n).*log(alpha_beta_exact(fi,n)./alpha_beta(fi,n)));   % DKL(Pmean(n),u_exact(n),Pvar(n),var_exact(n));
    f_alpha_beta_LT(n)= sum(alpha_beta_exact(fi_LT,n).*log(alpha_beta_exact(fi_LT,n)./alpha_beta_LT(fi_LT,n)));  % DKL(Pmean_LT(n),u_exact(n),Pvar_LT(n),var_exact(n));
    %     title(n);
    %     clf;
end

Pmean=PmeanT;
Pvar=PvarT;

Pmean_LT=Pmean_LTT;
Pvar_LT=Pvar_LTT;

BCJR.Pmean=Pmean;
BCJR.Pvar=Pvar;
BCJR.Taw=Taw;
BCJR.phase_shift=phase_shift;
% % BCJR.Pmean=u_exact.';
% % BCJR.Pvar=var_exact.';
Taw(pilot)=[];
Pmean(pilot)=[];
Pvar(pilot)=[];
Pmean_LT(pilot)=[];
Pvar_LT(pilot)=[];
alpha_gamma_beta_exact(pilot,:)=[];
y(pilot)=[];
y_LT=y.*exp(-1j*Pmean_LT);
y=y.*exp(-1j*Pmean);
variance_PN=Pvar;
variance_PN_LT=Pvar_LT;
XX=hMap.Xsort(:);
M=length(XX);

kk0=hMap.pntk0;
kk1=hMap.pntk1;

[K,kM]=size(hMap.pntk0);
LLR=zeros(length(y),kM);
LLR_LT=zeros(length(y),kM);
LLR_exact=zeros(length(y),kM);

L=zeros(length(y),kM);
L_LT=zeros(length(y),kM);
L_exact=zeros(length(y),kM);

xref    = repmat( XX.',length(y),1 );

yy      = repmat( y,1,M );
yy_LT      = repmat( y_LT,1,M );
Taw      = repmat( Taw,1,M );
vVar_PN = repmat( variance_PN, 1, M );
vVar_PN_LT = repmat( variance_PN_LT, 1, M );
%%% LT  %%%
y2 = yy; xref2 = xref;

LLR_log=-abs(y2-xref2).^2/N0;  %% Euclidean distance
LLR_log_LT=-abs(yy_LT-xref2).^2/N0;  %% Euclidean distance

%%% add correction

v2=N0+2*vVar_PN.*abs(xref2).^2;
v2_LT=N0+2*vVar_PN_LT.*abs(xref2).^2;

LLR_log=  LLR_log + imag(conj(xref2).*y2).^2*2.*(vVar_PN/N0)./v2 -0.5*log(v2);
LLR_log=LLR_log+ log(1+exp(log(Taw)-LLR_log));
LLR_log_LT=  LLR_log_LT + imag(conj(xref2).*y2).^2*2.*(vVar_PN_LT/N0)./v2_LT -0.5*log(v2_LT);
LLR_log_exact= log(alpha_gamma_beta_exact);
for k=1:kM
    LLR(:,k)=jac_log( LLR_log(:,kk1(:,k)) ) - jac_log( LLR_log(:,kk0(:,k)) );
    LLR_LT(:,k)=jac_log( LLR_log_LT(:,kk1(:,k)) ) - jac_log( LLR_log_LT(:,kk0(:,k)) );
    LLR_exact(:,k)=jac_log( LLR_log_exact(:,kk1(:,k)) ) - jac_log( LLR_log_exact(:,kk0(:,k)) );
end

% LR= exp( LLR_log );
% LR_LT=  exp( LLR_log_LT );
% LR_exact= alpha_gamma_beta_exact;
% for k=1:kM
%     L(:,k)=( sum(LR(:,kk1(:,k)),2) )./( sum(LR(:,kk0(:,k)),2) );
%     L_LT(:,k)=( sum(LR_LT(:,kk1(:,k)),2) ) ./ ( sum(LR_LT(:,kk0(:,k)),2) );
%     L_exact(:,k)=( sum(LR_exact(:,kk1(:,k)),2) ) ./ ( sum(LR_exact(:,kk0(:,k)),2) );
% end
%
% f_b=sum(L_exact.*log(L_exact./L));  %DKL(A,B,C,D);
% f_b_LT=sum(L_exact.*log(L_exact./L_LT));  %DKL(A,B,C,D);

P1=exp(LLR)./(1+exp(LLR));
P0=1./(1+exp(LLR));

P1_LT=exp(LLR_LT)./(1+exp(LLR_LT));
P0_LT=1./(1+exp(LLR_LT));

P1_exact=exp(LLR_exact)./(1+exp(LLR_exact));
P0_exact=1./(1+exp(LLR_exact));

f_1=sum(P1_exact.*log(P1_exact./P1));
f_1_LT=sum(P1_exact.*log(P1_exact./P1_LT));

f_0=sum(P0_exact.*log(P0_exact./P0));
f_0_LT=sum(P0_exact.*log(P0_exact./P0_LT));
fSP1=f_1+f_0;
fLT=f_1_LT+f_0_LT;

LLR_exact = LLR_exact';
LLR_exact = LLR_exact(:);
LLR_exact = -LLR_exact(1:Nc);
dec_llr_ph=step(hDec,LLR_exact);
dec_llr_ph=dec_llr_ph(1:Nb);
Binest_ph =  (dec_llr_ph<0) ;
difference     = ( Binest_ph ~= Data);
count_exact =(sum(difference));
PER_exact   = any( difference );

LLR = LLR';
LLR = LLR(:);
LLR = -LLR(1:Nc);
dec_llr_ph=step(hDec,LLR);
dec_llr_ph=dec_llr_ph(1:Nb);
Binest_ph =  (dec_llr_ph<0) ;
difference     = ( Binest_ph ~= Data);
count =(sum(difference));
PER   = any( difference );

LLR_LT = LLR_LT';
LLR_LT = LLR_LT(:);
LLR_LT = -LLR_LT(1:Nc);
dec_llr_ph=step(hDec,LLR_LT);
dec_llr_ph=dec_llr_ph(1:Nb);
Binest_ph =  (dec_llr_ph<0) ;
difference     = ( Binest_ph ~= Data);
count_LT =(sum(difference));
PER_LT   = any( difference );
[count_exact count count_LT]
end
