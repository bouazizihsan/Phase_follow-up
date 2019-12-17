function [BCJR]=Phase_track_BCJR_Gauss_Real_Phase_No_Class( y, N0, sigw2, X, pilot_symbols,pilot,Lframe,Phase_noise,is_BCJRP )

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
pdf_Gauss   = @(t,m,s2) exp( -abs(t-m ).^2./(2*s2) )./sqrt(2*pi*s2); % pdf(t|m), normpdf(x,mu,sigma.^2)
pdf_Gauss_log   = @(t,m,s2) ( -abs(t-m ).^2./(2*s2) )-0.5*log((2*pi*s2));
CC={'r--','b--','g--','k--'};
%processing
phase=angle(y(pilot,:)./pilot_symbols);
new_phases=unwrap(phase);
% %_____ au lieu de unwrap_________________________
% for n=1:length(phase)-1
%
%     deph= phase(n+1,:)-phase(n,:);
%     deph=deph(:);
%     fi=find(deph > pi); deph(fi)= deph(fi)-2.*pi;
%     fi=find(deph < -pi); deph(fi)= deph(fi)+2.*pi;
%     phase(n+1,:)=phase(n,:)+deph.';
% end
% %________________________________________________
%  % # Pilot

% newphases=zeros(size(new_phases));
% for n=1:length(new_phases)
% 
%     if(n==1)
%             sigman_1=Inf;
%             sigman1= (N0./(2*abs(pilot_symbols(n+1,:)).*abs(y(pilot(n)+Lframe,:))))+ Lframe.*sigw2;
%             sigman = (N0./(2*abs(pilot_symbols(n,:)).*abs(y(pilot(n),:)))) ;
%             newphases(n,:) =(new_phases(n,:)./sigman + new_phases(n+1,:)./sigman1)./(1./sigman+1./sigman_1+1./sigman1);
% 
%         elseif(n==length(pilot))
%             sigman1=Inf;
%             sigman_1=(N0./(2*abs(pilot_symbols(n-1,:)).*abs(y(pilot(n)-Lframe,:))))+ Lframe.*sigw2;
%             sigman = (N0./(2*abs(pilot_symbols(n,:)).*abs(y(pilot(n),:)))) ;
%             newphases(n,:) =(new_phases(n,:)./sigman + new_phases(n-1,:)./sigman_1)./(1./sigman+1./sigman_1+1./sigman1);
% 
%         else
%             sigman_1=(N0./(2*abs(pilot_symbols(n-1,:)).*abs(y(pilot(n)-Lframe,:))))+ Lframe.*sigw2;
%             sigman1= (N0./(2*abs(pilot_symbols(n+1,:)).*abs(y(pilot(n)+Lframe,:))))+ Lframe.*sigw2;
%             sigman = (N0./(2*abs(pilot_symbols(n,:)).*abs(y(pilot(n),:)))) ;
%             newphases(n,:)=(new_phases(n,:)./sigman + new_phases(n-1,:)./sigman_1 + new_phases(n+1,:)./sigman1)./(1./sigman+1./sigman_1+1./sigman1);
% 
%     end
% 
% 
% end
phase_shift=(interp1(pilot', new_phases, (1:Ns))).' ;
% phase_ref=Phase_noise;

%%%%%%%%%%%%%%%%%%   Interpolation Lineaire %%%%%%%%%%%%%%
%  new_phases=phase;
% phase_ref = zeros(Ns,1);
% var_ref = zeros(Ns,1);
% for n=1:length(pilot)
%     
% %     if(n==1)
% %         sigman_1=Inf;
% %         sigman1= (N0./(2*abs(pilot_symbols(n+1,:)).*abs(y(pilot(n)+Lframe,:))))+ Lframe.*sigw2;
% %         sigman = (N0./(2*abs(pilot_symbols(n,:)).*abs(y(pilot(n),:)))) ;
% % %         sigman=Inf;
% %         phase_ref(pilot(n),:) =(new_phases(n,:)./sigman + new_phases(n+1,:)./sigman1)./(1./sigman+1./sigman_1+1./sigman1);
% %         var_ref(pilot(n),1)=1./(1./sigman+1./sigman_1+1./sigman1);
% %     elseif(n==length(pilot))
% %         sigman1=Inf;
% %         sigman_1=(N0./(2*abs(pilot_symbols(n-1,:)).*abs(y(pilot(n)-Lframe,:))))+ Lframe.*sigw2;
% %         sigman = (N0./(2*abs(pilot_symbols(n,:)).*abs(y(pilot(n),:)))) ;
% % %         sigman=Inf;
% %         phase_ref(pilot(n),:) =(new_phases(n,:)./sigman + new_phases(n-1,:)./sigman_1)./(1./sigman+1./sigman_1+1./sigman1);
% %         var_ref(pilot(n),1)=1./(1./sigman+1./sigman_1+1./sigman1);
% %     else
% %         sigman_1=(N0./(2*abs(pilot_symbols(n-1,:)).*abs(y(pilot(n)-Lframe,:))))+ Lframe.*sigw2;
% %         sigman1= (N0./(2*abs(pilot_symbols(n+1,:)).*abs(y(pilot(n)+Lframe,:))))+ Lframe.*sigw2;
% %         sigman = (N0./(2*abs(pilot_symbols(n,:)).*abs(y(pilot(n),:)))) ;
% % %         sigman=Inf;
% %         phase_ref(pilot(n),:)=(new_phases(n,:)./sigman + new_phases(n-1,:)./sigman_1 + new_phases(n+1,:)./sigman1)./(1./sigman+1./sigman_1+1./sigman1);
% %         var_ref(pilot(n),1)=1./(1./sigman + 1./sigman_1 + 1./sigman1);
% %     end
%     
%     sigman = (N0./(2*abs(pilot_symbols(n,:)).^2)) ;
%        z=conj(y(pilot(n)))*exp(1j*new_phases(n));
%         phase_ref(pilot(n),:)=new_phases(n)-imag(z.*(pilot_symbols(n,:)))./abs(pilot_symbols(n,:)).^2;
%         var_ref(pilot(n),1)=sigman;
%     
% end
% 
% 
% for n=1:Ns % pL < n < (p+1)L
%     if( pilot_pos(n)==0)
%         ind=find(pilot<n);
%         nn= max(ind);
%         
%         %         (n-nn-(nn-1)*Lframe)
%         %         (nn+Lframe-n)
%         sigman_1=(N0./(2*abs(pilot_pos(pilot(nn),:)).*abs(y(pilot(nn),:))))+ (n-pilot(nn)).*sigw2;
%         sigman1= (N0./(2*abs(pilot_pos(pilot(nn+1),:)).*abs(y(pilot(nn+1),:))))+ (pilot(nn+1)-n).*sigw2;
%         
%         phase_ref(n,:)=(phase_ref(pilot(nn),:)./sigman_1 + phase_ref(pilot(nn+1),:)./sigman1)./(1./sigman_1+1./sigman1);
%         var_ref(n,1)=1./(1./sigman_1+1./sigman1);
%         
%     end
%     
% end
% BCJR.Pmean_apos1=phase_ref;
% BCJR.Pvar_apos1=var_ref;
% phase_ref(pilot)=[];
% var_ref(pilot)=[];
% BCJR.Pvar=var_ref;
% BCJR.Pmean=phase_ref;

%%%%%%%%%%% Fin d'Interpolation Lineaire %%%%%%%%%%%%%%
%processing
% phase=angle(y(pilot,:)./pilot_symbols);
% new_phases=unwrap(phase);
% phase_shift=(interp1(pilot', new_phases, (1:Ns))).' ;
y=y.*exp(-1j*phase_shift);
dd_phase_shift = diff( phase_shift ); 
phase_ref=(phase_shift-phase_shift);
for n=1:Ns
    z=conj(y(n))*exp(1j*phase_ref(n));
    if( pilot_pos(n)~=0)
        X_Hard=pilot_pos(n);
    else

        X_Hard=X;

    end
    Mnl_p0{1,n} =phase_ref(n)-imag(z.*(X_Hard))./abs(X_Hard).^2;
    Snl_Xci{1,n} = N0/2.*abs(X_Hard).^2;
    Snl{1,n} = N0./(2*abs(X_Hard).^2);
    Xci_p0 {1,n}= pdf_Gauss( real( z*X_Hard ) ,abs(X_Hard).^2, Snl_Xci{1,n});

    %     theta_bar=angle(y(n)./X_Hard)-Phase_noise(n);
    %
    %     Mnl_p0{1,n} =Phase_noise(n)-tan(theta_bar)./2;
    %     Snl{1,n} = N0./(2*abs(X_Hard.*y(n)).*cos(theta_bar));
    %     Xci_p0{1,n}= pdf_Gauss( y(n)  ,(X_Hard.*exp(1j*Phase_noise(n))), N0).*exp((sin(theta_bar).^2)./(2*cos(theta_bar)).*(2*abs(y(n).*X_Hard)./N0));
    %     gamma=(repmat(Xci_p0{1,n}.',Np,1)+pdf_Gauss_log(repmat(p,1,length(X_Hard)),repmat(Mnl_p0{1,n}.',Np,1),repmat(Snl{1,n}.',Np,1)));
    %     ggg=gamma(:);
    %     gamma=exp(gamma-jac_logH(ggg.'));
    %     figure(2)
    %     for i=1:length(X_Hard)
    %         hold on
    %         plot(p,gamma(:,i),CC{i});
    %
    %     end
    %     plot(Phase_noise(n),0,'r*');
    %     clf
end

%ALPHA
parameters1{1,1}=0;
parameters1{1,2}=inf;
parameters1{1,3}=log(1);
for n=2:Ns
    [Mnl,Sigma,ph]=Projection_of_product(Mnl_p0{1,n-1},Snl{1,n-1},log(Xci_p0{1,n-1}),0,parameters1{n-1,1},parameters1{n-1,2},parameters1{n-1,3},0,0,1,2,sigw2,Phase_noise(n));
    parameters1{n,1}=Mnl-dd_phase_shift(n-1);
    parameters1{n,2}=Sigma;
    parameters1{n,3}=ph;
end

%beta
parameters2{Ns,1}=0;
parameters2{Ns,2}=inf;
parameters2{Ns,3}=log(1);
for n=Ns-1:-1:1
    [Mnl,Sigma,ph]=Projection_of_product(Mnl_p0{1,n+1},Snl{1,n+1},log(Xci_p0{1,n+1}),0,parameters2{n+1,1},parameters2{n+1,2},parameters2{n+1,3},0,0,1,2,sigw2,Phase_noise(n));
    parameters2{n,1}=Mnl+dd_phase_shift(n);
    parameters2{n,2}=Sigma;
    parameters2{n,3}=ph;
end

Pmean=zeros(Ns,1);
Pvar=zeros(Ns,1);
for n=1:Ns
    if(n==Ns)
        [Pmean(n,:),Pvar(n,:),ph]=Projection_of_product(parameters1{n,1},parameters1{n,2},parameters1{n,3},0,parameters2{n,1},parameters2{n,2},parameters2{n,3},0,0,1,2,0,Phase_noise);
    else
        [Pmean(n,:),Pvar(n,:),ph]=Projection_of_product(parameters2{n,1},parameters2{n,2},parameters2{n,3},0,parameters1{n,1},parameters1{n,2},parameters1{n,3},0,0,1,2,0,Phase_noise);
    end
end

% Pmean(pilot)=[];
% Pvar(pilot)=[];

BCJR.Pmean=Pmean;
BCJR.Pvar=Pvar;
BCJR.Phase_Ref=phase_ref;
BCJR.phase_shift=phase_shift;

end
