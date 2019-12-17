function [BCJR]=Phase_track_Tikhonov_Hard1( y, N0, sigw2, X, pilot_symbols,pilot,Lframe,AlphaGammaProj )
pdf_Tikh=@(t,k,u) exp(k.*cos(t-u))./(2*pi.*besseli(0,k));
NClass=2;
D=1;
Ns           =length(y);
NP           =length(pilot);
pilot_pos(pilot)= pilot_symbols;
M            =length(X);
Npp           = 2^9;
% Phase_limit  = pi; %3*sqrt(sigw2*Lframe/4);  % we suppose the phase varies in the interval (-Phase_limit, Phase_limit)
% dpp          = 2.*Phase_limit./Np;
pp           = linspace(-pi,pi,Npp)';

% %processing
% phase=unwrap(angle(y(pilot)./pilot_symbols));
% phase_shift=interp1( pilot' , phase, (1:Ns)');
% dphase_shift=diff(phase_shift);
% y=y.*exp(-1j.*phase_shift);
[BCJRP]=PilotPhaseTrackTikhonov( y(pilot), N0, Lframe*sigw2 , pilot_symbols,Lframe); % Lframe*sigw2
z1{1,1}=0;
Z_alpha=zeros(Ns,D);
for n=2:Ns
    
    if(pilot_pos(n-1)~=0)
        
        X_Hard=pilot_pos(n-1);
        %         MM=1;
        %         z=(2/N0)*y(n-1)*conj(X_Hard); %% if N0-->0 abs(z)-->inf bessel(0,abs(z))=Inf
    else
        X_Hard=X;
        %         MM=M;
        %         z=repmat(z1{1,n-1},MM,1) + (2/N0)*y(n-1)*conj(X_Hard); %% if N0-->0 abs(z)-->inf bessel(0,abs(z))=Inf
    end
    MM=length(X_Hard);
    z=repmat(z1{1,n-1},MM,1) + (2/N0)*y(n-1)*conj(X_Hard); %% if N0-->0 abs(z)-->inf bessel(0,abs(z))=Inf
    lambda= log(besseli(0,abs(z),1))+abs(z)-(abs(X_Hard).^2)/N0;
    
    %     if(n==2)
    %         lambda=abs(z)-0.5*log(abs(z))-(abs(X_Hard).^2/N0);
    %     else
    %         lambda= (abs(z)-repmat(abs(z1{1,n-1}),MM,1))-0.5*log(abs(z)./abs(z1{1,n-1}))-(abs(X_Hard).^2)/N0;
    %     end
    lambda=exp(lambda-jac_logH(lambda.'));%% marginalisation
    
    
    switch AlphaGammaProj
        case 'OneSideCorr'
            if(n~=Ns)
                nn=min(find(pilot>n));
            else
                nn=NP;
            end
            %                 ZCorr=GaussToTikh(BCJRP.Z_beta(nn),(pilot(nn)-n).*sigw2); % (nn-1)*Lframe+1 = pilot(nn)
            ZCorr=gamma_fct((pilot(nn)-n).*sigw2,BCJRP.Z_beta(nn));
            %% Product Then Division
            z_gamma =gamma_fct(sigw2,z);
            [z_gamma,W]= TikhProduct(ZCorr,z_gamma);
            lambda=lambda.*exp(W);
            [zT] = CMVM(lambda,z_gamma,D);
            [zT,W]= TikhDivision(zT,ZCorr); %% Tikh(z_gamma)/Tikh(ZCorr)
            
            %             %% Mr. Leszek
            %             z_gamma =gamma_fct(sigw2,z);
            %             [zz,W]= TikhProduct(ZCorr,z_gamma);
            %             lambda_Prim=lambda.*exp(W);
            %             lambda=remove_irrelevant_terms( lambda,lambda_Prim);
            %             lambda=lambda./sum(lambda);
            %             [zT] = CMVM(lambda,z_gamma,D);
            
            %% Hsan
            %             if(MM>1)
            %                             d=TDKL(z,ZCorr) + TDKL(z,CMVM(lambda,z,D));
            %                             [val ind] = sort(d);
            %                             z=z(ind(1:NClass));
            %                             lambda=lambda(ind(1:NClass));
            %                             z_gamma =gamma_fct(sigw2,z);
            %                             [zT] = CMVM(lambda,z_gamma,D);
            %             else
            %                 zT=gamma_fct(sigw2,z);
            %             end
        case 'BothSideCorr'
            if(n~=Ns)
                nn=min(find(pilot>n));
            else
                nn=NP;
            end
            %                 Zbeta_pilot= GaussToTikh(BCJRP.Z_beta(nn),(pilot(nn)-n).*sigw2);
            Zbeta_pilot=gamma_fct((pilot(nn)-n).*sigw2,BCJRP.Z_beta(nn));
            %                 Zalpha_pilot=GaussToTikh(BCJRP.Z_alpha(nn-1),(n-pilot(nn-1)).*sigw2);
            Zalpha_pilot=gamma_fct((n-pilot(nn-1)).*sigw2,BCJRP.Z_alpha(nn-1));
            ZCorr= Zbeta_pilot+Zalpha_pilot;
            
            %%Mr. Leszek
            z_gamma =gamma_fct(sigw2,z);
            [zz,W]= TikhProduct(ZCorr,z_gamma);
            lambda_Prim=lambda.*exp(W);
            lambda=remove_irrelevant_terms( lambda,lambda_Prim);
            lambda=lambda./sum(lambda);
            
            [zT] = CMVM(lambda,z_gamma,D);
            %% Hsan
            
            %                 d=TDKL(z,ZCorr) + TDKL(z,CMVM(lambda,z,D));
            %                 [val ind] = sort(d);
            %                 z=z(ind(1:NClass));
            %                 lambda=lambda(ind(1:NClass));
            %                 z_gamma =gamma_fct(sigw2,z);
            %                 [zT] = CMVM(lambda,z_gamma,D);
            
        otherwise
            
            z_gamma =gamma_fct(sigw2,z);
            [zT] = CMVM(lambda,z_gamma,D);
    end
    
    if(any(isnan(zT)))
        stp=1;
    end
    z1{1,n}=zT;
    Z_alpha(n,1)=zT;
    
end

%Beta

z2{1,Ns}=0;
Z_beta=zeros(Ns,D);

for n=Ns-1:-1:1
    if(pilot_pos(n+1) ~= 0)
        X_Hard=pilot_pos(n+1);
    else
        X_Hard=X;
    end
    MM=length(X_Hard);
    z=repmat(z2{1,n+1},MM,1) + (2/N0)*y(n+1)*conj(X_Hard);
    lambda= log(besseli(0,abs(z),1))+abs(z)-(abs(X_Hard).^2)/N0;
    %     if(n==(Ns-1))
    %         lambda=abs(z)-0.5*log(abs(z))-(abs(X_Hard).^2/N0);
    %     else
    %         lambda=(abs(z)-repmat(abs(z2{1,n+1}),MM,1))-0.5*log(abs(z)./abs(z2{1,n+1}))-(abs(X_Hard).^2)/N0;
    %     end
    lambda=exp(lambda-jac_logH(lambda.'));
    
    switch AlphaGammaProj
        case 'OneSideCorr'
            if(n~=1)
                nn=max(find(pilot<n));
            else
                nn=1;
            end
            %                 ZCorr=GaussToTikh(BCJRP.Z_alpha(nn),(n-pilot(nn)).*sigw2);  % (nn-1)*Lframe+1 = pilot(nn)
            ZCorr=gamma_fct((n-pilot(nn)).*sigw2,BCJRP.Z_alpha(nn));
            %% Product Then Division
%             z_gamma =gamma_fct(sigw2,z);
%             [z_gamma,W]= TikhProduct(ZCorr,z_gamma);
%             lambda=lambda.*exp(W);
%             [zT] = CMVM(lambda,z_gamma,D);
%             [zT,W]= TikhDivision(zT,ZCorr); %% Tikh(z_gamma)/Tikh(ZCorr)
            %% Mr. Leszek
                        z_gamma =gamma_fct(sigw2,z);
                        [zz,W]= TikhProduct(ZCorr,z_gamma);
                        lambda_Prim=lambda.*exp(W);
                        lambda=remove_irrelevant_terms( lambda,lambda_Prim);
                        lambda=lambda./sum(lambda);
                        [zT] = CMVM(lambda,z_gamma,D);
            %%Hsan
            %             if(MM>1)
            %                             d=TDKL(z,ZCorr) + TDKL(z,CMVM(lambda,z,D));
            %                             [val ind] = sort(d);
            %                             z=z(ind(1:NClass));
            %                             lambda=lambda(ind(1:NClass));
            %                             z_gamma =gamma_fct(sigw2,z);
            %                             [zT] = CMVM(lambda,z_gamma,D);
            %             else
            %                 zT=gamma_fct(sigw2,z);
            %             end
        case 'BothSideCorr'
            if(n~=1)
                nn=max(find(pilot<n));
            else
                nn=1;
            end
            %                 Zalpha_pilot=GaussToTikh(BCJRP.Z_alpha(nn),(n-pilot(nn)).*sigw2);
            Zalpha_pilot=gamma_fct((n-pilot(nn)).*sigw2,BCJRP.Z_alpha(nn));
            %                 Zbeta_pilot=GaussToTikh(BCJRP.Z_beta(nn),(pilot(nn+1)-n).*sigw2);
            Zbeta_pilot=gamma_fct((pilot(nn+1)-n).*sigw2,BCJRP.Z_beta(nn+1));
            ZCorr= Zbeta_pilot+Zalpha_pilot;
            %% Mr. Leszek
            z_gamma =gamma_fct(sigw2,z);
            [zz,W]= TikhProduct(ZCorr,z_gamma);
            lambda_Prim=lambda.*exp(W);
            lambda=remove_irrelevant_terms( lambda,lambda_Prim);
            lambda=lambda./sum(lambda);
            
            [zT] = CMVM(lambda,z_gamma,D);
            %%Hsan
            %
            %                 d=TDKL(z,ZCorr) + TDKL(z,CMVM(lambda,z,D));
            %                 [val ind] = sort(d);
            %                 z=z(ind(1:NClass));
            %                 lambda=lambda(ind(1:NClass));
            
        otherwise
            z_gamma =gamma_fct(sigw2,z);
            [zT] = CMVM(lambda,z_gamma,D);
    end
    if(any(isnan(zT)))
        stp=1;
    end
    z2{1,n}=zT;
    Z_beta(n,1)=zT;
    
    
end

BCJR.Z_beta=Z_beta; %.*exp(1j.*phase_shift);
BCJR.Z_alpha=Z_alpha;% .*exp(1j.*phase_shift);

end

% function z=gamma_fct(sig,z)
% 
% z=z./(1+abs(z).*sig);
% 
% end

function z=gamma_fct(sig,z)
k1=1/sig;
k2=abs(z);
u=angle(z);
circular_var2=1-(besseli(1,k1,1)./besseli(0,k1,1) .* besseli(1,k2,1)./besseli(0,k2,1) );
k=1./(2*circular_var2);

z=k.*exp(1j*u);

end

function d=TDKL(A,B)
AA=abs(A);
BB=abs(B);
ua=angle(A);
ub=angle(B);
d=log(besseli(0,AA,1)./besseli(0,BB,1))+(AA)-(BB)+ besseli(1,BB,1)./besseli(0,BB,1).*(BB - AA.*cos(ua-ub));
%     d=AA.*(1-cos(ub-ua))-0.5.*log(AA./BB)+AA./(2*BB).*cos(ub-ua)  ;
end
function z=GaussToTikh(A,sigw2)
AA=abs(A);
ua=angle(A);
z=(1./((1./AA)+sigw2)).*exp(1j.*ua);

end
function [z,W]=TikhProduct(Z1,Z2)
% Z1= Zcorr
z= Z1 + Z2;
W= log(besseli(0,z,1))+abs(z)-log(besseli(0,Z2,1))-abs(Z2); % Weights

end
function [z,W]=TikhDivision(Z1,Z2)

z= Z1 - Z2;
W= abs(z)-abs(Z1)+abs(Z2)-0.5.*log(abs(z).*abs(Z2)./(abs(Z1))); % Weights

end
function phi_tilde=remove_irrelevant_terms( phi_tilde,phi_tilde_tmp)

CUT_threshold = 0.95; %% used if AlphaBetaPilotFilter='Y' to remove unnecessary elements from gamma
%PICK_f=0;

%% remove the irrelevant terms

[vv, ii] = sort( phi_tilde_tmp/sum(phi_tilde_tmp) , 'descend');
vv_cum = cumsum(vv);
[fi]=find( vv_cum > CUT_threshold );
ii_zero = ii( fi(2:end) );
phi_tilde( ii_zero ) = 0;
end