function [BCJR]=Phase_track_BCJR_Gauss3( y, N0, sigw2, hMap, pilot_pos, pilot_symbols, LLR_apr, PAR)

%%% phase tracking using BCJR
% alpha, beta are modelled as Gaussian; gamma as a SUM of Gaussians
% y         : received signal
% N0        : variance of the additive Gaussian noise
% sigw2     : variance of the phase noise
% X         : signal constellation
% pilot_pos : pointers to the pilot positions
% pilot_symbols : pilot symbols
% LLR_apr    : a priori LLR for all the bit (meaningful for the payload)
% PAR.Gamma_approx: which phasse reference is used for Gaussian approximation for the gamma elements
%               'LTmax' -use max and Taylor around max,
%               'LTmax-F' -use max and Taylor around max but reject peaks far from zero
%               'LT0' - use Taylor around zero, 'BLT0' -  use Taylor+bilinear around zero
%               'BLT0-F' % use biblinear transform around max but reject peaks far from zero
% PAR.AlphaBetaReboot: 'Y', restart phase estimation at each pilot
% PAR.AlphaBetaPilotFilter: 'Y'; use future/past pilots to filter out contributions from alpha/beta
%

% Last modification: July 2019

Ns  = length(y);
X = hMap.X;
if nargin<7, LLR_apr=[]; end
if nargin<8,
    GammaApproximation='SP';
    AlphaBetaPilotFilter='None';
else
    GammaApproximation=PAR.GammaApproximation;
    AlphaBetaPilotFilter=PAR.AlphaBetaPilotFilter;
end

PRIOR = ~isempty(LLR_apr);
if PRIOR
    Lbin = hMap.Lbin;
    mm = size(hMap.pntk0,2);
end
Npilots = length(pilot_symbols);
Lframe   = (Ns-1)/(Npilots-1); % duration of the frame (payload+ one pilot)
is_pilot = zeros(Ns,1);
is_pilot(pilot_pos) = pilot_symbols; %% pilot value at the pilot position
nn = rem( (0:Ns-1)'  ,Lframe ) ;  %% position relative the the previous pilot
previous_pilot = ceil( (1:Ns)'/Lframe ); %% index of the previous pilot
previous_pilot(end) = previous_pilot(end-1);

pdf_ypx     = @(y,p,x) exp( -abs(y-x.*exp( 1j*p) ).^2/N0 ); % pdf(y|p,x)
pdf_Gauss   = @(y,m,s2) exp( -(y-m ).^2./(2*s2) )./sqrt(2*pi*s2); % pdf(y|p,x) real Gaussian
diff_phase  = @(x) x + 2*pi.*(  (x< -pi) - (x>pi) ) ;
%pdf_Gauss   = @(y,m,s2) exp( -diff_phase(y-m).^2./(2*s2) )./sqrt(2*pi*s2); % pdf(y|p,x) real Gaussian

%%%%% REMOVE THIS PART =>
% define phases used to find the contributing symbols
%Phase_limit  =  3*sqrt(sigw2*Lframe/4 + N0/2);  %% we suppose the phase varies in the interval (-Phase_limit, Phase_limit)
%dpp =   Phase_limit/100;%hMap.phi_min;
%pp  =   (-Phase_limit:dpp:Phase_limit)';
%%%%% <= REMOVE THIS PART

%% Preprocessing: de-rotate the received signal
%Phase = angle(y(pilot_pos)./pilot_symbols); % crude phase from pilots

BCJR_pilots = Phase_from_pilots( y( pilot_pos ) , pilot_symbols, N0, Lframe*sigw2 );  %% filtered phase from pilots
Phase = BCJR_pilots.apo(:,1);
%Phase = BCJR_pilots.gamma(:,1);

phase_shift = interp1( pilot_pos, Phase, (1:Ns)');
dd_phase_shift = diff( phase_shift ); %% difference in the argument of the pdf

y  = y.*exp( -1j* phase_shift );

%%% partial aposteriori variance of the pilot using all other past(future) pilots
if 1
    var_pilot_beta = 1./( 1./BCJR_pilots.beta(:,2)+ 1./BCJR_pilots.gamma(:,2) );  % estimate of the variance of the future pilot
    mean_pilot_beta = ( BCJR_pilots.beta(:,1)./BCJR_pilots.beta(:,2)...
        + BCJR_pilots.gamma(:,1)./BCJR_pilots.gamma(:,2) ).*var_pilot_beta;
    var_pilot_alpha = 1./( 1./BCJR_pilots.alpha(:,2)+ 1./BCJR_pilots.gamma(:,2) );  % estimate of the variance of the future pilot
    mean_pilot_alpha = ( BCJR_pilots.alpha(:,1)./BCJR_pilots.alpha(:,2)...
        + BCJR_pilots.gamma(:,1)./BCJR_pilots.gamma(:,2) ).*var_pilot_alpha;
    
    
    v_alpha_n = var_pilot_alpha( previous_pilot ) + nn*sigw2;
    v_beta_n = var_pilot_beta( previous_pilot+1 )+ ( Lframe-nn )*sigw2;
                   
    var_pilot_n = 1./( 1./v_alpha_n + 1./v_beta_n );
                    
    m_alpha_n = mean_pilot_alpha( previous_pilot );
    m_beta_n = mean_pilot_beta( previous_pilot+1 );
    mean_pilot_n = ( m_alpha_n./v_alpha_n + m_beta_n./v_beta_n ) .* var_pilot_n;

else % just use next pilot
    var_pilot_alpha = BCJR_pilots.gamma(nnpilot,2);
    mean_pilot_alpha = BCJR_pilots.gamma(nnpilot,1); %%
    var_pilot_beta  = var_pilot_alpha;
    mean_pilot_beta = mean_pilot_alpha;
end

%% parameters
alpha = zeros(Ns,2); % first column contains means, second - variance
beta  = zeros(Ns,2);
gamma = zeros(Ns,2);
pdf_W = [0, sigw2];

%% ----------G A M M A
for n=1:Ns
    if is_pilot(n)~=0
        xhard = is_pilot( n );
    else
        xhard      =  X;  %%% all constellation symbols
    end
    switch GammaApproximation
        case 'SP' % use Laplace approximation around max of the phase
            g_m{1,n}   = angle( y(n)./xhard );
            g_s2{1,n}  = N0./( 2*abs(xhard.*y(n)) );
            g_xi{1,n}  = pdf_Gauss( abs( y(n)), abs(xhard ), N0/2 )./ sqrt(abs(xhard.*y(n))) ;
        case 'LTmax' % use Taylor series around max of the phase
            g_m{1,n}   = angle( y(n)./xhard );
            g_s2{1,n}  = N0./( 2*abs(xhard).^2 );
            g_xi{1,n}  = pdf_Gauss( abs( y(n)), abs(xhard ), N0/2 )./ abs(xhard) ;
        case 'LT' % use Taylor series around phase=zero
            g_m{1,n}   = imag( y(n)./xhard );
            g_s2{1,n}  = N0./( 2*abs(xhard).^2 );
            if real(is_pilot(n))~=0  %% change calculation for pilots
               g_m{1,n}   = angle( y(n)./xhard ); 
               g_s2{1,n}  = N0./( 2*abs(xhard.*y(n)) );
            end
            g_xi{1,n}  = exp( -abs(y(n)-xhard).^2/N0 + imag( y(n)./xhard).^2 .*abs(xhard).^2 /N0 ) ./ abs(xhard) ;
    end
    
    if real(is_pilot(n))==0 && PRIOR
        Proba_xhard=  exp( Lbin( dd(ia,n), 1:mm )*(LLR_apr(n,:)') );
        g_xi{1,n} = g_xi{1,n} .*Proba_xhard;
    end
    
end



%% A L P H A alpha(:,n)
alpha(1,:) = [0 inf]; % no past for the first symbol; inf= non-informative prior, ie uniform distribution
for n=2:Ns-1
    if n==2 %% after pilot positions
        phi_tilde = 1;
        m_tilde   = g_m{n-1};
        s2_tilde   = g_s2{n-1};
    else
        %parametrize the product \alpha_{n-1} * \gamma_{n-1}
        phi_tilde = g_xi{n-1} .* pdf_Gauss( alpha(n-1,1), g_m{n-1}, g_s2{n-1} + alpha(n-1,2)  );
        m_tilde   =  ( alpha(n-1,1)*g_s2{n-1} + g_m{n-1}*alpha(n-1,2) )./( alpha(n-1,2) + g_s2{n-1} );
        s2_tilde  =  1./( 1./alpha(n-1,2) + 1./g_s2{n-1} );
        
        %%% OPTIONAL FILTERING 
        % compare alpha with prior knowledge from future pilots and remove irrelevant terms
        if real(is_pilot(n))==0 % only at payload positions
            switch AlphaBetaPilotFilter
                case 'OneSideCut'
                    var_pilot_nn = v_beta_n( n ); %% variance of the estimate from the future pilot
                    mean_pilot_nn = m_beta_n( n ) - phase_shift(n) ;   %% mean of the estimate from the future pilot
                    [phi_tilde] = remove_irrelevant_terms( m_tilde, phi_tilde, mean_pilot_nn, var_pilot_nn );%+ s2_tilde );
                case 'BothSidesCut'
                    var_pilot_nn = var_pilot_n( n );
                    mean_pilot_nn = mean_pilot_n( n )- phase_shift(n);
                    [phi_tilde] = remove_irrelevant_terms( m_tilde, phi_tilde, mean_pilot_nn, var_pilot_nn );
                case 'Multiply'
                    var_pilot_nn = v_beta_n( n ); %% variance of the estimate from the future pilot
                    mean_pilot_nn = m_beta_n( n ) - phase_shift(n) ;   %% mean of the estimate from the future pilot
                    phi_tilde = phi_tilde .* pdf_Gauss( m_tilde, mean_pilot_nn, var_pilot_nn + s2_tilde );
                case 'EP'
                    
            end
        end
    end
    
    % project the product on the Gaussian
    sum_phi =  sum(phi_tilde);
    mu_check  =  ( m_tilde' * phi_tilde )/sum_phi;
    sig2_check=  ( ( s2_tilde + m_tilde.^2)' * phi_tilde )/sum_phi - mu_check^2;
    
    alpha(n,:) = [mu_check , sig2_check] + [ -dd_phase_shift(n-1) , sigw2 ]; %+ pdf_W ;
end


%% B E T A beta(:,n)
beta(Ns,:) = [0, inf]; % no future for the last symbol
for n=Ns-1:-1:2
    if  n==(Ns-1)
        phi_tilde = 1;
        m_tilde   = g_m{n+1};
        s2_tilde   = g_s2{n+1};
    else
        phi_tilde = g_xi{n+1} .* pdf_Gauss( beta(n+1,1), g_m{n+1}, g_s2{n+1} + beta(n+1,2)  );
        m_tilde   =  ( beta(n+1,1)*g_s2{n+1} + g_m{n+1}*beta(n+1,2) )./( beta(n+1,2) + g_s2{n+1} );
        s2_tilde  =  1./( 1./beta(n+1,2) + 1./g_s2{n+1} );
        
        %%% OPTIONAL FILTERING 
        if real(is_pilot(n))==0 % only at payload positions
            switch AlphaBetaPilotFilter
                case 'OneSideCut'
                    var_pilot_nn = v_alpha_n( n ); %% variance of the estimate from the future pilot
                    mean_pilot_nn = m_alpha_n( n )- phase_shift(n) ;   %% mean of the estimate from the future pilot
                    [phi_tilde] = remove_irrelevant_terms( m_tilde, phi_tilde, mean_pilot_nn, var_pilot_nn );%+ s2_tilde );
                case 'BothSidesCut'
                    var_pilot_nn = var_pilot_n( n );
                    mean_pilot_nn = mean_pilot_n( n )- phase_shift(n);
                    [phi_tilde] = remove_irrelevant_terms( m_tilde, phi_tilde, mean_pilot_nn, var_pilot_nn );
                case 'Multiply'
                    var_pilot_nn = v_beta_n( n ); %% variance of the estimate from the future pilot
                    mean_pilot_nn = m_beta_n( n ) - phase_shift(n) ;   %% mean of the estimate from the future pilot
                    phi_tilde = phi_tilde .* pdf_Gauss( m_tilde, mean_pilot_nn, var_pilot_nn + s2_tilde );

            end
        end
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
BCJR.var  = TT(:,2);
BCJR.alpha = alpha;
BCJR.beta  = beta;
BCJR.gamma = gamma;

%BCJR.pp    = pp;
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
[~,fi]=min(abs(x-m));
y(fi)=1;

function phi_tilde=remove_irrelevant_terms( m_tilde, phi_tilde, mean_pilot_nn, var_pilot_nn )

diff_phase  = @(x) x + 2*pi.*(  (x< -pi) - (x>pi) ) ;
pdf_Gauss   = @(y,m,s2) exp( -diff_phase(y-m).^2./(2*s2) )./sqrt(2*pi*s2); % pdf(y|p,x) real Gaussian


CUT_threshold = 0.95; %% used if AlphaBetaPilotFilter='Y' to remove unnecessary elements from gamma
%PICK_f=0;

%% remove the irrelevant terms
phi_tilde_tmp = phi_tilde .* pdf_Gauss( m_tilde, mean_pilot_nn, var_pilot_nn  );

[vv, ii] = sort( phi_tilde_tmp/sum(phi_tilde_tmp) , 'descend');
vv_cum = cumsum(vv);
[fi]=find( vv_cum > CUT_threshold );
ii_zero = ii( fi(2:end) );
phi_tilde( ii_zero ) = 0;

%             switch PICK_f
%                 case 0%% the terms contributing to the distribution survive
%                     [fi]=find( vv_cum > CUT_threshold );
%                     ii_zero = ii( fi(2:end) );
%                 case 1%% only the maximum term survives
%                     ii_zero = ii( 2:end );
%                 case 2%%  the terms with the smallest mean survive
%                     [vv, ii]=sort( abs(m_tilde-mean_pilot_nn), 'ascend' );
%                     ii_zero = ii( 3:end );
%             end

