function [BCJR]=Phase_track_BCJR_Gauss_D( y, N0, sigw2, X, pilot_symbols,pilot,Lframe,payload,Gamma_Apprx,Pilot_Gamma_Apprx,AlphaGammaProj)
BCJRP=[];
ss           =size(y);
NClass=2;
D            =ss(2);
Ns           =ss(1);
pilot_pos    = zeros(Ns,2);
pilot_pos(pilot,:)= pilot_symbols;
M            =length(X);
N0_2           =N0.*eye(D);
sig2        =sigw2.*eye(D);

[BCJRP]=BCJR_PilotD( y(pilot,:), N0, Lframe*sig2, pilot_symbols,Lframe,payload,Pilot_Gamma_Apprx);

% Np           = 2^8;
% Npilots      = length(pilot_symbols);
% Phase_limit  =3*sqrt(sig2*Lframe/4);  % we suppose the phase varies in the interval (-Phase_limit, Phase_limit)
%
% dpp         = 2.*pi./Np;
% p           = linspace(-pi,pi,Np).';

% [M1,M2]=meshgrid(p);
% % M1=M1.';
% % M2=M2.';
% M11=M1(:);
% M22=M2(:);
% pp=[M11 M22];

pdf_Gauss=@(t,m,s2,d) real(exp(-0.5*sum(((t-m )'*inv(s2).*(t-m ).'),2) )./((2*pi).^(d/2)*det(s2).^(0.5)));

pdf_Gauss_log=@(t,m,s2,d) ((-0.5*(t-m )'*inv(s2)*(t-m )) - log( (2*pi).^(d/2)*det(s2).^(0.5) ));

% g= (reshape(pdf_Gauss( pp.' ,0 ,sigw2,D),Np,Np));

%processing
new_phases=BCJRP.apo; % unwrap(angle(y(pilot,:)./pilot_symbols)); % BCJRP.mean_apos;%
%     newphases=zeros(size(new_phases));
%     for n=1:length(new_phases)
%
%         if(n==1)
%             sigman_1=Inf;
%             sigman1= (N0./(2*abs(pilot_symbols(n+1,:)).*abs(y(pilot(n)+Lframe,:))))+ Lframe.*sig2;
%             sigman = (N0./(2*abs(pilot_symbols(n,:)).*abs(y(pilot(n),:)))) ;
%             newphases(n,:) =(new_phases(n,:)./sigman + new_phases(n+1,:)./sigman1)./(1./sigman+1./sigman_1+1./sigman1);
%
%         elseif(n==length(pilot))
%             sigman1=Inf;
%             sigman_1=(N0./(2*abs(pilot_symbols(n-1,:)).*abs(y(pilot(n)-Lframe,:))))+ Lframe.*sig2;
%             sigman = (N0./(2*abs(pilot_symbols(n,:)).*abs(y(pilot(n),:)))) ;
%             newphases(n,:) =(new_phases(n,:)./sigman + new_phases(n-1,:)./sigman_1)./(1./sigman+1./sigman_1+1./sigman1);
%
%         else
%             sigman_1=(N0./(2*abs(pilot_symbols(n-1,:)).*abs(y(pilot(n)-Lframe,:))))+ Lframe.*sig2;
%             sigman1= (N0./(2*abs(pilot_symbols(n+1,:)).*abs(y(pilot(n)+Lframe,:))))+ Lframe.*sig2;
%             sigman = (N0./(2*abs(pilot_symbols(n,:)).*abs(y(pilot(n),:)))) ;
%             newphases(n,:)=(new_phases(n,:)./sigman + new_phases(n-1,:)./sigman_1 + new_phases(n+1,:)./sigman1)./(1./sigman+1./sigman_1+1./sigman1);
%
%         end
%
%
%     end

phase_shift=[(interp1(pilot', new_phases(:,1), (1:Ns)')) interp1(pilot', new_phases(:,2), (1:Ns)')];
%phase_shift=0*phase_shift;
dd_phase_shift = diff( phase_shift );
% y=y.*exp(-1j*phase_shift);
phase_ref=phase_shift; %(phase_shift-phase_shift);%

Theta_zero=zeros(M,D,Ns);
for n=1:Ns
    auxi=[];
    aux_SNL=[];
    aux1=[];
    %     z=conj(y(n,:)).*exp(1j*phase_ref(n,:));
    if( norm(pilot_pos(n,:))~=0)
        X_Hard=pilot_pos(n,:);
    else
        
        X_Hard=X;
        
    end
    MM=size(X_Hard,1);
    switch Gamma_Apprx
        case 'SP'
%             Moy{1,n} =angle(repmat(y(n,:),MM,1)./X_Hard);
            m=[unwrap([repmat(phase_shift(n,1),1,MM);(angle(y(n,1)./(X_Hard(:,1)))).']) unwrap([repmat(phase_shift(n,2),1,MM);(angle(y(n,2)./(X_Hard(:,2)))).'])];
            Moy{1,n}=[m(2,1:MM).' m(2,MM+1:end).'];
            for i=1:MM
                % |X|*|Y|
                aux1=[aux1;diag(N0./(2*abs(X_Hard(i,:)).^2) )];
                aux_SNL=[aux_SNL;diag(N0./(2*abs(X_Hard(i,:)).*abs(y(n,:))))] ;
                auxi=[auxi;(log(1./prod(( abs(X_Hard(i,:)).*(sqrt(abs(X_Hard(i,:).*y(n,:)))) ))) + pdf_Gauss_log( abs(y(n,:)./(X_Hard(i,:))).',1,aux1((i-1)*D+1:i*D,:),D) )];
            end
            
        case 'LTmax'
%             Moy{1,n} =angle(repmat(y(n,:),MM,1)./X_Hard);
            m=[unwrap([repmat(phase_shift(n,1),1,MM);(angle(y(n,1)./(X_Hard(:,1)))).']) unwrap([repmat(phase_shift(n,2),1,MM);(angle(y(n,2)./(X_Hard(:,2)))).'])];
            Moy{1,n}=[m(2,1:MM).' m(2,MM+1:end).'];
            for i=1:MM
                % |X|*|Y|
                aux1=[aux1;diag(N0./(2*abs(X_Hard(i,:)).^2) )];
                aux_SNL=[aux_SNL;diag(N0./(2*abs(X_Hard(i,:)).^2))] ;
                auxi=[auxi;(log(1./prod(( abs(X_Hard(i,:)).^2))) + pdf_Gauss_log( abs(y(n,:)./(X_Hard(i,:))).',1,aux1((i-1)*D+1:i*D,:),D) )];
            end
        case 'LT'
            z=conj(y(n,:)).*exp(1j*phase_ref(n,:));
            Moy{1,n} =repmat(phase_ref(n,:),MM,1)-imag(repmat(z,MM,1).*(X_Hard))./abs(X_Hard).^2;
            for i=1:MM
                % |X|^2
                aux_SNL=[aux_SNL;diag(N0./(2*abs(X_Hard(i,:)).^2) )];
                auxi=[auxi;(pdf_Gauss_log( (real(z.*X_Hard(i,:))).',(abs(X_Hard(i,:)).^2).',aux_SNL((i-1)*D+1:i*D,:),D) )];
                
            end
        otherwise
            error('unknown data Gamma Apprx method ');
            
    end
    
    Theta_zero(:,:,n)=[Moy{1,n};zeros(M-size(Moy{1,n},1),ss(2))];
    Xci_p0{1,n}=auxi;
    Snl{1,n}=aux_SNL;
    
end


%ALPHA
% alpha=zeros(Np,Np,Ns);
% alpha(:,:,1)=1;
alpha_p0{1,1}=[1 1] ;
alpha_p0{1,2}=diag(ones(1,D).*Inf);
alpha_p0{1,3}=log(1); % les poids sont en log.
% corr=1;
alpha_p0{Ns,1}=0;
alpha_p0{Ns,2}=diag(ones(1,D));
alpha_p0{Ns,3}=0;
for n=2:Ns-1
    switch AlphaGammaProj
        case 'OneSideCorr'
            
            nn=min(find(pilot>n));
            m=unwrap([phase_shift(n,:);BCJRP.parameters2{1,nn}]);
            Corr{1,1}=m(2,:)-phase_shift(n,:);
            Corr{2,1}=BCJRP.parameters2{2,nn}+(pilot(nn)-n).*sigw2;
            
        case 'BothSideCorr'
            
            nn=min(find(pilot>n));
            m=unwrap([phase_shift(n,:);BCJRP.parameters2{1,nn}]);
            beta_pilot{1,1}=m(2,:)-phase_shift(n,:);
            beta_pilot{2,1}= BCJRP.parameters2{2,nn}+(pilot(nn)-n).*sigw2;
            
            m=unwrap([phase_shift(n,:);BCJRP.parameters1{1,nn-1}]);
            alpha_pilot{1,1}=m(2,:)-phase_shift(n,:);
            alpha_pilot{2,1}=BCJRP.parameters1{2,nn-1}+(n-pilot(nn-1)).*sigw2;
            [Corr{1,1},Corr{2,1},Corr{3,1}]=Produit_Alpha_Gamma_D(alpha_pilot{1,1},alpha_pilot{2,1},0,beta_pilot{1,1},beta_pilot{2,1},0,D);
        otherwise
            Corr=0;
    end
    [cor,Mnl,S,Poids]=Projection_of_productD(Moy{1,n-1},Snl{1,n-1},Xci_p0{1,n-1},alpha_p0{n-1,1},alpha_p0{n-1,2},alpha_p0{n-1,3},Corr,AlphaGammaProj,D);
    alpha_p0{n,1}= Mnl-repmat(dd_phase_shift(n-1,:),size(Mnl,1),1);
    alpha_p0{n,2}=S+repmat(sig2,size(Poids,1),1);
    alpha_p0{n,3}= Poids;
    
end

%Beta
% beta=zeros(Np,Np,Ns);
% beta(:,:,Ns)=1;
beta_p0{Ns,1}=[1 1] ;
beta_p0{Ns,2}=diag(ones(1,D).*Inf);
beta_p0{Ns,3}=0;
% corr=1;
beta_p0{1,1}=0;
beta_p0{1,2}=diag(ones(1,D));
beta_p0{1,3}=0;

for n=Ns-1:-1:2
    switch AlphaGammaProj
        case 'OneSideCorr'
            nn=max(find(pilot<n));
            m=unwrap([phase_shift(n,:);BCJRP.parameters1{1,nn}]);
            Corr{1,1}=m(2,:)-phase_shift(n,:);
            Corr{2,1}=BCJRP.parameters1{2,nn}+(n-pilot(nn)).*sigw2;
            
        case 'BothSideCorr'
            nn=max(find(pilot<n));
            m=unwrap([phase_shift(n,:);BCJRP.parameters1{1,nn}]);
            alpha_pilot=[m(2,:)-phase_shift(n,:);BCJRP.parameters1{2,nn}+(n-pilot(nn)).*sigw2];
            
            m=unwrap([phase_shift(n,:);BCJRP.parameters2{1,nn+1}]);
            beta_pilot=[m(2,:)-phase_shift(n,:);BCJRP.parameters2{2,nn+1}+(pilot(nn+1)-n).*sigw2];
            [Corr{1,1},Corr{2,1},Corr{3,1}]=Produit_Alpha_Gamma_D(alpha_pilot(1,:),alpha_pilot(2,:),0,beta_pilot(1,:),beta_pilot(2,:),0,D);
        otherwise
            Corr=0;
    end
    [cor,Mnl,S,Poids]=Projection_of_productD(Moy{1,n+1},Snl{1,n+1},Xci_p0{1,n+1},beta_p0{n+1,1},beta_p0{n+1,2},beta_p0{n+1,3},Corr,AlphaGammaProj,D);
    beta_p0{n,1}= Mnl+repmat(dd_phase_shift(n,:),size(Mnl,1),1);
    beta_p0{n,2}=S+repmat(sigw2,size(Poids,1),1);
    beta_p0{n,3}= Poids;
    
    
end

m_alpha_beta=zeros(Ns,D);
for i=1:Ns
    if(i==Ns)
        [ind,Mnl,S,ph]=Projection_of_productD(alpha_p0{i,1},alpha_p0{i,2},alpha_p0{i,3},beta_p0{i,1},beta_p0{i,2},beta_p0{i,3},0,'None',D);
        
    else % cas ou il y a pas de varience inf
        [ind,Mnl,S,ph]=Projection_of_productD(beta_p0{i,1},beta_p0{i,2},beta_p0{i,3},alpha_p0{i,1},alpha_p0{i,2},alpha_p0{i,3},0,'None',D);
        
    end
    
    alpha_beta{i,1}= Mnl+phase_shift(i,:);
    alpha_beta{i,2}= S;
    alpha_beta{i,3}=ph;
    %     alpha_beta{i,1}= [Mnl;zeros(NClass^2-size(Mnl,1),D)];
    %     alpha_beta{i,2}= [S;zeros((NClass^2-size(Mnl,1))*2,D)];
    %     alpha_beta{i,3}=[ph;-Inf.*ones(NClass^2-size(ph,1),1)];
    
    
end


BCJR.uapox=[alpha_beta{:,1}];
m=[alpha_beta{payload,1}].';
v=[alpha_beta{payload,2}].';
p=[alpha_beta{payload,3}].';
alpha_beta=[];
for i=1:length(p)
    alpha_beta{i,1}=m((i-1)*D+1:i*D,:).';
    alpha_beta{i,2}=v((i-1)*D+1:i*D,:).';
    alpha_beta{i,3}=p(i,:).';
end
BCJR.P_Ref=Theta_zero(:,:,payload);%alpha_beta{:,1};%phase_ref(payload,:);%
BCJR.alpha_beta=alpha_beta;

end


% %Alpha
% alpha_p1{1,1}=[0 0] ;
% alpha_p1{1,2}=diag(ones(1,D).*Inf);
% corr=1;
% for n=2:Ns
%     if( (norm(pilot_pos(n-1,:))~=0) && ((1-corr)>10^(-4)))
%         aux1=[alpha_p1{n-1,1};(0*(1:D))];
%         aux2=[alpha_p1{n-1,2};diag(Inf*(1:D))];
%         aux3=[log(corr);log(1-corr)];
%         corr=1;
%         [ind,Mnl,S,PHI]=Projection_of_productD(Moy{1,n-1},Snl{1,n-1},Xci_p0{1,n-1},aux1,aux2,aux3,0,1,D,M1,M2,Phase_noise(n,:));
%
%     else
%         [ind,Mnl,S,PHI]=Projection_of_productD(Moy{1,n-1},Snl{1,n-1},Xci_p0{1,n-1},alpha_p1{n-1,1},alpha_p1{n-1,2},0,0,1,D,M1,M2,Phase_noise(n,:));
%     end
%     %     [ind,Mnl,S,PHII]=Projection_of_productD(Moy{1,n-1},Snl{1,n-1},Xci_p0{1,n-1},alpha_p0{n-1,1},alpha_p0{n-1,2},0,D,M1,M2,Phase_noise(n,:));
%     %     [ind,Mnll,SS,PHI]=Projection_of_productD(Moy{1,n-1},Snl{1,n-1},Xci_p0{1,n-1},beta_p0{n-1,1},beta_p0{n-1,2},0,D,M1,M2,Phase_noise(n,:));
%     % [exp(PHII-jac_logH(PHII.')) exp(PHI-jac_logH(PHI.'))]
%
%     [cor,Mnll,SS,PHI]=Projection_of_productD(Mnl,S,PHI,beta_p0{n-1,1},beta_p0{n-1,2},0,1,1,D,M1,M2,Phase_noise(n,:));
%     [ind,Mnll,SS,PHI]=Projection_of_productD(beta_p0{n-1,1},beta_p0{n-1,2},0,Mnll,SS,0,1,1,D,M1,M2,Phase_noise(n,:));
%
%     % %      [ind,Mnlll,SSS,PHI]=Projection_of_productD(Mnl,S,0,beta_p0{n-1,1},beta_p0{n-1,2},2,1,D,M1,M2,Phase_noise(n,:));
%
%     %       [alpha_p0{n,1};Mnl]
%     %           Mnl=Mnl(ind,:);
%     %     S = S((ind-1)*D+1:ind*D,:);
%     if((det(SS)>0).* ((diag(SS).')>[0 0]))
%         alpha_p1{n,1}= Mnll;
%         alpha_p1{n,2}=SS+sigw2;
%         corr=corr*cor;
%     else
%         alpha_p1{n,1}= alpha_p0{n,1};
%         alpha_p1{n,2}=alpha_p0{n,2};
%     end
%
%
%
%     %     alpha(:,:,n)= ifft2( fft2( alpha(:,:,n-1).* gamma(:,:,n-1), Np,Np ) .* fft2(fftshift(g), Np,Np)  );
%     %     alpha(:,:,n)=alpha(:,:,n)./(sum(sum(alpha(:,:,n))));
%     %     figure(9)
%     %     contour(M1,M2,(alpha(:,:,n)));
%     %     hold on
%     %     plot(Phase_noise(n,1),Phase_noise(n,2),'r*');
%     %     hold off
%     %     figure(10)
%     %     contour(M1,M2,gamma(:,:,n-1).*alpha(:,:,n-1));
%     %     hold on
%     %     plot(Phase_noise(n,1),Phase_noise(n,2),'r*');
%     %     hold off
% end
%
% %Beta
% beta_p1{Ns,1}=[0 0] ;
% beta_p1{Ns,2}=diag(ones(1,D).*Inf);
% corr=1;
% for n=Ns-1:-1:1
%     if( (norm(pilot_pos(n+1,:))~=0) && ((1-corr)>10^(-4)))
%         aux1=[beta_p1{n+1,1};(0*(1:D))];
%         aux2=[beta_p1{n+1,2};diag(Inf*(1:D))];
%         aux3=[log(corr);log(1-corr)];
%         corr=1;
%         [ind,Mnl,S,PHI]=Projection_of_productD(Moy{1,n+1},Snl{1,n+1},Xci_p0{1,n+1},aux1,aux2,aux3,0,1,D,M1,M2,Phase_noise(n,:));
%
%     else
%     %     [ind,Mnl,S,PHII]=Projection_of_productD(Moy{1,n+1},Snl{1,n+1},Xci_p0{1,n+1},beta_p0{n+1,1},beta_p0{n+1,2},0,D,M1,M2,Phase_noise(n,:));
%         [ind,Mnl,S,PHI]=Projection_of_productD(Moy{1,n+1},Snl{1,n+1},Xci_p0{1,n+1},beta_p1{n+1,1},beta_p1{n+1,2},0,0,1,D,M1,M2,Phase_noise(n,:));
%     %      [ind,Mnll,SS,PHI]=Projection_of_productD(Moy{1,n+1},Snl{1,n+1},Xci_p0{1,n+1},alpha_p0{n+1,1},alpha_p0{n+1,2},0,D,M1,M2,Phase_noise(n,:));
%     end
%     %     Mnl=Mnl(ind,:);
%     %     S = S((ind-1)*D+1:ind*D,:);
%
%     [cor,Mnll,SS,PHIII]=Projection_of_productD(Mnl,S,PHI,alpha_p0{n+1,1},alpha_p0{n+1,2},0,2,1,D,M1,M2,Phase_noise(n,:));
%     [ind,Mnll,SS,PHIII]=Projection_of_productD(alpha_p0{n+1,1},alpha_p0{n+1,2},0,Mnll,SS,0,1,2,D,M1,M2,Phase_noise(n,:));
%
%     % [exp(PHII-jac_logH(PHII.')) exp(PHI-jac_logH(PHI.')) exp(PHIII-jac_logH(PHIII.'))]
%     %     Mnl=Mnl(ind,:);
%     %     S = S((ind-1)*D+1:ind*D,:);
%     if( (det(SS)>0).* ((diag(SS).')>[0 0]))
%         if(size(SS)~=[2 2])
%             stp=1;
%         end
%         beta_p1{n,1}= Mnll;
%         beta_p1{n,2}=SS+sigw2;
%         corr=corr*cor;
%     else
%         beta_p1{n,1}= beta_p0{n,1};
%         beta_p1{n,2}=beta_p0{n,2};
%     end
%
%     %      beta(:,:,n)= ifft2( fft2( beta(:,:,n+1).* gamma(:,:,n+1), Np,Np ) .* fft2(fftshift(g), Np,Np)  );
%     %      beta(:,:,n)=beta(:,:,n)./(sum(sum(beta(:,:,n))));
%     %     figure(9)
%     %     contour(M1,M2,(beta(:,:,n)));
%     %     hold on
%     %     plot(Phase_noise(n,1),Phase_noise(n,2),'r*');
%     %     hold off
%     %
%     %     figure(10)
%     %     contour(M1,M2,gamma(:,:,n+1).*beta(:,:,n+1));
%     %     hold on
%     %     plot(Phase_noise(n,1),Phase_noise(n,2),'r*');
%     %     hold off
%     %     fig=(pdf_Gauss([M11.';M22.'] ,repmat((beta_p0{n+1,1}).',1,length(M11)),beta_p0{n+1,2},D ));
%     %     fig=reshape(fig,[],length(M1));
%     %     fig=fig./(sum(sum(fig)));
%     %     figure(11)
%     %     contour(M1,M2,gamma(:,:,n+1).*fig);
%     %     hold on
%     %     plot(Phase_noise(n,1),Phase_noise(n,2),'r*');
%     %     hold off
% end

% p_apo=repmat(beta.*alpha,1,1,1,M).*gammaX;
% p_apo=p_apo./repmat(sum(sum(sum(p_apo,1),2),4),Np,Np,1,M);
% p_apo=sum(sum(p_apo));

% % for i=1:Ns
% %     Pmean(i,2)=angle(exp(1j.*p).'*sum(p_apo(:,:,i),2).*dpp);
% %     Pmean(i,1)=angle(exp(1j.*p).'*sum(p_apo(:,:,i),1)'.*dpp);
% % end

%Product of Alpha and Beta ==> one gaussian
% beta_alpha=zeros(Np,Np,Ns);