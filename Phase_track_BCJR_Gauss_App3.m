function [BCJR]=Phase_track_BCJR_Gauss_App3( y, N0, sigw2, hMap, pilot_pos, pilot_symbols)

%%% phase tracking using BCJR
% alpha, beta are modelled as Gaussian; gamma as a SUM of Gaussians
% y         : received signal
% N0        : variance of the additive Gaussian noise
% sigw2     : variance of the phase noise
% X         : signal constellation
% pilot_pos : pointers to the pilot positions
% pilot_symbols : pilot symbols
  
Ns  = length(y);
X = hMap.X;
Npilots = length(pilot_symbols);
Lframe   = (Ns-1)/(Npilots-1); % duration of the frame (payload+ one pilot)
is_pilot = zeros(Ns,1);
is_pilot(pilot_pos) = pilot_symbols; %% pilot value at the pilot position
Np=2^9;
% define phases used to find the contributing symbols
Phase_limit  =  3*sqrt(sigw2*Lframe/4);  %% we suppose the phase varies in the interval (-Phase_limit, Phase_limit)
dpp =  2.*Phase_limit./Np; ;%hMap.phi_min;
pp  =  linspace(-Phase_limit,Phase_limit,Np)';
%Np  =   length(pp);

%% Preprocessing: derotate the received signal
Phase = angle(y(pilot_pos)./pilot_symbols); % phase from pilots
%--
% unwrapping the phase 
dPhase= diff(Phase); %% phase difference; suppose it's in (-pi, pi)
fi=find( dPhase>pi ); dPhase( fi ) = dPhase( fi ) -2*pi ;
fi=find( dPhase<-pi ); dPhase( fi ) = dPhase( fi ) +2*pi ;
Phase_lin = Phase(1) + [0;cumsum( dPhase)];
phase_shift = interp1( pilot_pos, Phase_lin, (1:Ns)');
dd_phase_shift = diff( phase_shift ); %% difference in the argument of the pdf

y   =  y.*exp( -1j* phase_shift );

%% parameters
alpha = zeros(Ns,2); % first column contains means, second - variance
beta  = zeros(Ns,2);
gamma = zeros(Ns,2);
pdf_W = [0, sigw2];

pdf_ypx     = @(y,p,x) exp( -abs(y-x.*exp( 1j*p) ).^2/N0 ); % pdf(y|p,x)
pdf_Gauss   = @(y,m,s2) exp( -abs(y-m ).^2./(2*s2) )./sqrt(2*pi*s2); % pdf(y|p,x) real Gaussian

%% ----------G A M M A
YY = repmat(y.',Np,1); PP = repmat(pp,1,Ns);
VV = YY.*exp( -1i*PP );
[XX_hard, dd] = hard_QAM( VV(:), X );

%%LS===== used for testing.... actual values of gamma( p )
%XX_hard = reshape( step( hMod, dd ), Np, Ns );
XX_hard = reshape( XX_hard, Np, Ns );

%LS======

dd      = reshape( dd, Np, Ns);  %% indices to the symbols X_hard

for n=1:Ns
    if is_pilot(n)~=0
        xhard = is_pilot( n );
        %g_m{1,n}   = imag( y(n)./xhard ); % mean
        g_m{1,n}   = angle( y(n)./xhard ); % mean
        g_s2{1,n}  = N0./(2*abs(xhard).^2); % variance
        g_xi{1,n}  = pdf_Gauss( abs( y(n)./xhard ) ,1, g_s2{1,n} )./ abs(xhard).^2 ;
    else
        [ddu{1,n},ia] = unique(dd(:,n)); %% indices to unique symbols x_hard
        %xhard    = step( hMod, ddu{1,n});
        if (any(ia ~= 1))
            s=1;
        end
        xhard      =  X( dd(ia,n) );
        g_m{1,n}   = angle( y(n)./xhard );
        g_s2{1,n}  = N0./(2*abs(xhard).^2);
        g_xi{1,n}  = pdf_Gauss( abs( y(n)./xhard ) ,1, g_s2{1,n} )./ abs(xhard).^2 ;
  
    end
%     Lgg=length(g_m{n});
%     semilogy(pp,gammaL(:,n),'b',pp,...
%         pdf_Gauss( repmat(pp,1,Lgg), repmat(g_m{n}',Np,1), repmat(g_s2{n}',Np,1) ) * g_xi{n},'r'),grid
end




%% A L P H A alpha(:,n) <-  alpha_{n-1}(p_n)
alpha(1,:) = [0 inf]; % no past for the first symbol; inf= non-informative prior, ie uniform distribution
alphaL(1,:) = [0 inf];
for n=2:Ns
    if n==2, 
        phi_tilde = 1;
        m_tilde   = g_m{n-1};
        s2_tilde   = g_s2{n-1};
    else
        %parametrize the product \alpha_{n-1} * \gamma_{n-1}
        phi_tilde = g_xi{n-1} .* pdf_Gauss( alpha(n-1,1), g_m{n-1}, g_s2{n-1} + alpha(n-1,2)  );
        m_tilde   =  ( alpha(n-1,1)*g_s2{n-1} + g_m{n-1}*alpha(n-1,2) )./( alpha(n-1,2) + g_s2{n-1} );
        s2_tilde  =  1./( 1./alpha(n-1,2) + 1./g_s2{n-1} );
    end
    % project the product on the Gaussian
    sum_phi =  sum(phi_tilde);
    mu_check  =  ( m_tilde' * phi_tilde )/sum_phi;
    sig2_check=  ( ( s2_tilde + m_tilde.^2)' * phi_tilde )/sum_phi - mu_check^2;
    
    alpha(n,:) = [mu_check , sig2_check] + [ -dd_phase_shift(n-1) , sigw2 ]; %+ pdf_W ;
end


%% B E T A beta(:,n) <- beta_{n+1}(p_n)
beta(Ns,:) = [0, inf]; % no future for the last symbol
betaL(Ns,:) = [0, inf];
for n=Ns-1:-1:1
    if n==(Ns-1),
        phi_tilde = 1;
        m_tilde   = g_m{n+1};
        s2_tilde   = g_s2{n+1};
    else
        phi_tilde = g_xi{n+1} .* pdf_Gauss( beta(n+1,1), g_m{n+1}, g_s2{n+1} + beta(n+1,2)  );
        m_tilde   =  ( beta(n+1,1)*g_s2{n+1} + g_m{n+1}*beta(n+1,2) )./( beta(n+1,2) + g_s2{n+1} );
        s2_tilde  =  1./( 1./beta(n+1,2) + 1./g_s2{n+1} );
    end
    % project the product on the Gaussian
    sum_phi =  sum( phi_tilde );
    mu_check  =  ( m_tilde' * phi_tilde )/sum_phi;
    sig2_check=  ( ( s2_tilde + m_tilde.^2)' * phi_tilde )/sum_phi - mu_check^2;

    beta(n,:) = [mu_check, sig2_check] + [ dd_phase_shift(n) , sigw2 ]; %+ pdf_W ;
end

%% Combine
TT = mult_Gauss( alpha, beta );

BCJR.Pmean = TT(:,1);
BCJR.Pvar  = TT(:,2);
BCJR.alpha = alpha;
BCJR.beta  = beta;
BCJR.gamma = gamma;

BCJR.pp    = pp;
BCJR.phase_shift = phase_shift;

function c = mult_Gauss( a, b )
% parameters of the Gaussian after multiplication
for n=1:size(a,1)
if isfinite(a(n,2)) && isfinite(b(n,2)) %% finite variance
    ss = a(n,2) + b(n,2);
    c(n,1) = ( a(n,1).*b(n,2) +  b(n,1).*a(n,2) )./ss;
    c(n,2) = ( b(n,2).*a(n,2) )./ss;
else
    if isinf(a(n,2))
        c(n,:) = b(n,:); % a is non-informative
    else
        c(n,:) = a(n,:); % b is non-informative
    end
end
end

function c = get_Gauss( pp, pdf )
% obtain the estimates of the mean and the variance
dpp=pp(2)-pp(1);
pdf=pdf/sum(pdf)/dpp;
c(1,1) = ( pp'*pdf )'*dpp; % estimated mean
c(1,2) = ( (pp.^2)'*pdf )'*dpp - (c (1) ).^2; % estimated variance

function y=pdf_Dirac(x, m)
% discrete Dirac impulse
y=zeros(length(x),1);
[zz,fi]=min(abs(x-m));
y(fi)=1;
