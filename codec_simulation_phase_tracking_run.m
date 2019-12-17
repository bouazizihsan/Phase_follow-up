function  [ BER, PER, GMI, PHI_stat ]= codec_simulation_phase_tracking_run( hChan, hEnc, hDec, SNRdB )

%%% perform the whole simulation: coding, transmission, phase tracking and decoding

%% setup the parameters
Npilots = hChan.Npilots; %% number of pilots
Lframe  = hChan.Lframe;  %% payload + one pilot
Lpayload= Lframe-1;       % number of symbols between the pilots
Ns      = (Npilots-1)*(Lframe) +1;      % number of symbols
pilots  = (1: Lframe: Ns)';         % positions of the pilots (assume uniformly distributed)
payload = (1:Ns)';
payload(pilots) = []; %% indices to the payload ONLY

sigw2   = hChan.sigw2;    % variance of the phase noise
%SNRdB   = hChan.SNRdB;    % SNR in dB
SNRlin  = 10.^(SNRdB/10); % Linear SNR
N0      = 1./SNRlin;

M       = hChan.M; % M-QAM
m       = log2(M);
hMap    = hChan.hMap;
[K,mm]  = size(hMap.pntk0);

[hMapPilot] = create_constellation(2,'Gray','QAM'); % pilots constellation is 4-QAM
modSignal   = zeros(Ns,1);

%%% coding info
[Nk,Nc]=size(hEnc.ParityCheckMatrix);
Nb=Nc-Nk; %% number of information bits

LLR2GMI = @(x) log2( 1+ exp( -abs(x) ) ) + x.*(x>0)*log2(exp(1)) ;


BER = 0;
PER = 0;
PHI_stat = zeros(Ns,2);
GMI = zeros(length(payload),1);

Nblocks = hChan.Nblocks; % number of independent blocks
Nblocks_stop = hChan.Nblocks_stop; 
for jj=1:Nblocks
    if mod(jj,ceil(Nblocks/10))==0, fprintf('.'); end  %% prints out when with 10% steps
    %% Generate the data
    % random binary data
    data  = logical( randi([0 1], Nb, 1) );
    % Encode LDPC
    encodedData    = step(hEnc, data);
    % pad with random bits to fit multiple of mm
    encodedDataPad  = [encodedData;randi([0 1], rem(Nc,mm), 1)];
    % interleaving
    % INT = randperm(Nc);
    % encodedDataPad = encodedDataPad(INT);
    
    % modulate
    Ncodedsymbols = length(encodedDataPad)/mm; %% number of symbols bearing the coded data = (Ns-Npilots)
    WholeData = reshape( encodedDataPad, mm, Ncodedsymbols); %% reshape
    if mm~=m  %% if the number of coded bits is smaller than the number of bits in the constellation (Ungerboeck)
        WholeData = [WholeData;  logical(randi([0 1], m-mm, Ncodedsymbols))]; %% add random bits (uncoded)
    end
    encodeDataDec = bi2de( WholeData', 'left-msb' ); %% make it decimal
    modSignal(payload)  = hMap.Xsort(encodeDataDec+1);
    
    %% Generate the pilots
    Xpilots     = hMapPilot.X( randi( [1 4], Npilots, 1 ) ); % random 4-QAM pilot symbols
    modSignal(pilots) = Xpilots; % symbols
    
    %% transmission over the channel
    w   =   randn(Ns,1);    % iid phase noise
    Ww  =   sqrt(sigw2)*cumsum(w);     % Wiener process
    Z   =   sqrt(N0/2)*( randn(Ns,1) + 1i*randn(Ns,1) );    % Generate Ns complex random variables with  variance (complex) N0
    y   =   modSignal .* exp( 1i*Ww ) + Z;
    
    %% Tracking + signal rotation
    switch hChan.PhaseEstimation
        case 'P-Direct-Linear'
            [BCJR] = phase_rotation( y, N0, sigw2, pilots, Xpilots, 0); %0: direct reading from pilots
            y_shift = y.*exp( -1j* (BCJR.phase_shift) );
        case 'Artificial-Parabola' %% iid phase noise with the same marginal distribution (variance is a parabola) as after linear interpolation
            [BCJR] = phase_rotation( y, N0, sigw2, pilots, Xpilots, 0); 
            Wwparabola = w.*sqrt( BCJR.var ); %% iid phase noise with the same distribution as the original one
            y_shift = modSignal .* exp( 1i*Wwparabola ) + Z;
        case 'P-Direct-Linear-maxvar'  %%% uses the maximum phase-noise variance estimate at all positions
            [BCJR] = phase_rotation( y, N0, sigw2, pilots, Xpilots, 0); %0: direct reading from pilots
            variance_PN_max = sigw2*Lframe/4 + N0/4;
            BCJR.var = ones(Ns,1)*variance_PN_max;
            y_shift = y.*exp( -1j* (BCJR.phase_shift) );
        case 'P-BCJR-Linear' %% pilot estimates using BCJR, the rest is done via linear interpolation
            [BCJR] = phase_rotation( y, N0, sigw2, pilots, Xpilots, 1); %1: BCJR on pilots
            y_shift = y.*exp( -1j* (BCJR.phase_shift) );
        case 'P-BCJR-EM' %% pilot estimates using BCJR, the rest is done via linear interpolation + Expectation maximization
            [BCJR] = phase_rotation_EM( y, N0, sigw2, pilots, Xpilots, hMap);
            y_shift = y.*exp( -1j* (BCJR.phase_shift) );
        case 'All-Pilots' %% all symbols are used to estimate the phase 
            [BCJR] = phase_rotation( y, N0, sigw2, (1:Ns), modSignal, 2); % performance limit
            y_shift = y.*exp( -1j* (BCJR.phase_shift) );
        case 'P-Perfect-Linear' %% no additive noise at the pilots= phase is perfectly estimated
            y_temp = modSignal .* exp( 1i*Ww );
            [BCJR] = phase_rotation( y_temp, 0, sigw2, pilots, Xpilots, 0); %pilot phase perfectly known
            y_shift = y.*exp( -1j* (BCJR.phase_shift) );            
        case 'All-BCJR' % BCJR for all symbols,
            [BCJR] = Phase_track_BCJR_Gauss3( y, N0, sigw2, hMap, pilots, Xpilots, [], hChan.BCJR_PAR);
            %[BCJR] = Phase_track_BCJR_Gauss4( y, N0, sigw2, hMap, pilots, Xpilots);
            y_shift = y.*exp( -1j* (BCJR.phase_shift + BCJR.Pmean) );
    end
    
    %% phase variance
%     PHI_stat(:,1) = PHI_stat(:,1) + (BCJR.phase_shift - Ww);
%     PHI_stat(:,2) = PHI_stat(:,2) + (BCJR.phase_shift - Ww).^2;
    
    %% LLR calculations for the payload
    LLR   = LLR_PN( y_shift(payload), 1, N0, BCJR.var(payload), hMap, hChan.LLR_METHOD );
    
    %%%% MI for the indicated bits
    LLR_flip=LLR.*( 1-2*WholeData(1:mm,:)');
    % calculate the GMI per symbol
    GMI = GMI + sum( LLR2GMI( LLR_flip ), 2 )/mm;
    
    LLR   = LLR';
    LLR   = - LLR(:);%% negative sign to accomodate matlab   
    
    % deinterleaving
    % LLR(INT) = LLR;
    
    % Decoding
    dec_out_LLR   = step(hDec, LLR);
  
% Error Count on info bits
if 1
    receivedBits   = ( dec_out_LLR(1:Nb)<0 ); %% hard decision on info bits and remove padding
    difference     = ( receivedBits ~= data); 
    BER        = BER + sum( difference )/Nb;
else% Error Count on code bits
    %receivedBits   = ( dec_out_LLR(1:Nc)<0 ); %% hard decision on code bits
    receivedBits   = ( (dec_out_LLR(1:Nc)-LLR)<0 ); %% hard decision on extrinsic LLRs
    difference     = ( receivedBits ~= encodedData); 
    BER        = BER + sum( difference )/Nc;
end
    
    PER        = PER + any( difference );

    % find optimal parameter "s" for GMI
%     LLR_flip=LLR.*( 2*encodedData -1 );
%     ss  =  (0.95:0.01:1.15);
%     [GMImin,si] = min( sum(LLR2GMI( LLR_flip*ss ),1) ,[],2 );
%     ss(si)

if PER>=Nblocks_stop, break, fprintf('early stop\n'), end
end
Nblocks = jj;  %% take early stop into account
GMI     =   1 - GMI/ Nblocks;
BER     =   BER/Nblocks;
PER     =   PER/Nblocks;
PHI_stat=   PHI_stat/Nblocks;



