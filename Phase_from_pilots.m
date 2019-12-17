function [BCJR]=Phase_from_pilots( y_ref, x_ref, N0, sigw2)

%%% phase tracking using BCJR assume all the symbols are known
% y_ref     : received signal
% x_ref     : transmitted symbols
% N0        : variance of the additive Gaussian noise
% sigw2     : variance of the phase noise (from one sample to another)
%
% alpha, beta, gamma are modelled each as a Gaussian

Ns = length(y_ref);

%% Preprocessing: derotate the received signal
Phase_ref = angle( y_ref./x_ref ); % phase from symbols
% unwrapping the phase
Phase_ref = unwrap(Phase_ref);

gamma = [Phase_ref, N0./(2*abs(x_ref).^2)]; % [mean, variance]
%gamma = [Phase_ref, N0./(2*abs(x_ref.*y_ref))]; % [mean, variance] %% modification done 5-07-2019
alpha = zeros(Ns,2); % first column contains means, second - variance
beta  = zeros(Ns,2);

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


function c = mult_Gauss( a, b )
% parameters of the Gaussian after multiplication
%for n=1:size(a,1)
    c(:,2) = 1./( 1./b(:,2)+ 1./a(:,2) );
    c(:,1) = ( a(:,1)./a(:,2) +  b(:,1)./b(:,2) ).*c(:,2);
%end
