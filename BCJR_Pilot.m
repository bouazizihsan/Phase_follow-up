function [BCJR]=BCJR_Pilot( Y, N0, sigw2, X,Pilot_Gamma_Apprx)
Ns=length(Y);
switch Pilot_Gamma_Apprx
    case 'SP'
        gamma = [unwrap(angle(Y./X)) N0./(2*abs(X).*abs(Y))];
%         Xci_p0=pdf_Gauss( abs(Y) ,abs(X_Hard),N0/2)./ (sqrt(abs(X_Hard*Y)) );
    case 'LTmax'     
        gamma = [unwrap(angle(Y./X))  N0./(2*abs(X).^2)];
%         Xci_p0=pdf_Gauss( abs(Y) ,abs(X_Hard),N0/2)./(abs(X_Hard)) ;
end
%% A L P H A alpha(:,n) <-  alpha_{n-1}(p_n)
alpha(1,:) = [0 inf]; % no past for the first symbol; inf= non-informative prior, ie uniform distribution
for n=2:Ns
    alpha(n,:) = mult_Gauss( alpha(n-1,:) , gamma(n-1,:) ) + [0, sigw2];
end

%% B E T A beta(:,n) <- beta_{n+1}(p_n)
beta(Ns,:) = [0, inf]; % no future for the last symbol
for n=Ns-1:-1:1
    beta(n,:) = mult_Gauss( beta(n+1,:) , gamma(n+1,:) ) + [0, sigw2];
end

%% Combine
TT_ext = mult_Gauss( alpha, beta );
TT_apo = mult_Gauss( TT_ext, gamma);
BCJR.alpha = alpha;
BCJR.beta  = beta;
BCJR.gamma = gamma; % direct reading of the phase
BCJR.ext = TT_ext;  %% extrinsic (to get the estimate excluding the pilot)
BCJR.apo = TT_apo;  %% aposteriori (to get the estimate for the pilots)
BCJR.beta_gamma=mult_Gauss( gamma, beta );
BCJR.alpha_gamma=mult_Gauss( gamma, alpha );

function c = mult_Gauss( a, b )
% parameters of the Gaussian after multiplication
%for n=1:size(a,1)
    c(:,2) = 1./( 1./b(:,2)+ 1./a(:,2) );
    c(:,1) = ( a(:,1)./a(:,2) +  b(:,1)./b(:,2) ).*c(:,2);
%end
end
end
