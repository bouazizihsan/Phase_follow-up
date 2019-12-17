function [count,PER]=codec_simulation_PN_run( hChan, hEnc, hDec )


% perform the whole simulation: coding, transmission, phase tracking and decoding
% LS, March 2017, modified March 2018

%% setup the parameters

sigw2   = hChan.sigw2;    % variance of the phase noise
SNRdB   = hChan.SNRdB;    % SNR in dB
SNRlin  = 10.^(SNRdB/10); % Linear SNR
N0      = 1./SNRlin;
BCJR2   =[];
BCJR    =[];
M       = hChan.M; % M-QAM
m       = log2(M);
hMap    = hChan.hMap;
XCons=hMap.X(:);
[K,mm]  = size(hMap.pntk0);
PhaseEstimation=hChan.PhaseEstimation;
Lframe= hChan.Lframe;
%%% coding info
[Nk,Nc]=size(hEnc.ParityCheckMatrix);
Nb=Nc-Nk; %% number of information bits

%%% modulation info
m=log2(length(hMap.X));
[K,mm]=size(hMap.pntk0); %% mm number of bits per symbol
if K~=length(hMap.X)/2, error('wrong size of XOUT.pntk0'); end

LLR2GMI = @(x) log2( 1+ exp( -abs(x) ) ) + x.*(x>0)*log2(exp(1)) ; % avoid calculation of exp() with large arguments

count  = 0; %%% average bit error rate
count_tkh  = 0;
BER=zeros(2,1);
PER  = 0; %%% average bit packeteerror rate
BlocksN = hChan.BlocksN; %% number of transmitted blocks
SUM_METHOD=hChan.SUM_METHOD;
LLR_METHOD=hChan.LLR_METHOD;
%STAT = zeros(BlocksN,1);  %% register the number of bits in error in each block
%GMI  = zeros(BlocksN,1);  %% calculate the GMI at the decoder's output
for jj=1:BlocksN
    if mod(jj,ceil(BlocksN/10))==0, fprintf('.'); end  %% prints out when with 10% steps
    % random data
    data  = logical(randi([0 1], Nb, 1));
    % Encode LDPC
    encodedData    = step(hEnc, data);
    % pad with random bits to fit multiple of mm
    encodedDataPad  = [encodedData;randi([0 1], mm-rem(Nc,mm), 1)];
    % interleave
    
    % Modulate
    Ncodedsymbols=length(encodedDataPad)/mm; % number of symbols
    WholeData = reshape( encodedDataPad, mm, Ncodedsymbols); %% reshape
    if mm~=m  %% if the number of coded bits is smaller than the number of bits in the constellation (Ungerboeck)
        WholeData = [WholeData;  logical(randi([0 1], m-mm, Ncodedsymbols))]; %% add random bits (uncoded)
    end
    encodeDataDec = bi2de( WholeData', 'left-msb' ); %% make it decimal
    modSignal      = hMap.Xsort(encodeDataDec+1);
    
    pilots  = (1:Lframe: length(modSignal))';
    pilots=[pilots;length(modSignal)];
    Npilots=length( pilots );
    % Generate the pilots
    var=randi([1 M], Npilots, 1);
%     bits=de2bi( var, 'left-msb' );
    Xpilots     = hMap.X(var); % random 4-QAM pilot symbols
    modSignal(pilots) = Xpilots; % symbols
    % Channel
    Z  =  sqrt(N0/2)*(randn(Ncodedsymbols,1)+1i*randn(Ncodedsymbols,1));    % Generate ns complex random variables with  variance (complex) N0
    
    switch PhaseEstimation
        case 'yes'
            w   =   sqrt(sigw2).*randn(size(modSignal))+ 0; % sigma_p.^2: variance et 0= mean
            Phase_noise  =   cumsum(w);     % Wiener process
            receivedSignal  =   modSignal.*exp(1i*Phase_noise)+Z;
           [BCJR]  =Phase_track_Tikhonov( receivedSignal, N0, sigw2, XCons.', Xpilots,pilots,Lframe );
           %[BCJR2]    = Phase_track_BCJR_Gauss_AppL( receivedSignal, N0, sigw2, XCons.', Xpilots,pilots,Lframe );
        
        
        case 'None'
            Phase_noise=sqrt(sigw2)*randn(size(modSignal));
            receivedSignal  =   modSignal.*exp(1i*Phase_noise)+Z;
        otherwise
            error('unknown PhaseEstimation');
    end
    
    
    % Demodulate
    %LLR   =   LLR_PN(receivedSignal, 1, N0, sigw2, hMap, hChan.LLR_METHOD,SUM_METHOD);

    LLR_tikh =LLR_Tikhonov(M,receivedSignal(payload),XCons,N0,BCJR.Z_alpha,BCJR.Z_beta,BCJR.Fact_alpha,BCJR.Fact_beta,SUM);
      
    LLR     =   LLR_cluster_Gauss(M,receivedSignal(payload),XCons.',N0,sigw2,BCJR2.L_alpha,BCJR2.L_beta,BCJR2.parameters1,BCJR2.parameters2,SUM_METHOD,LLR_METHOD);
        
    LLR   =   LLR';
    LLR   =   - LLR(:);%% negative sign to accomodate matlab
    LLR   =   LLR( 1: Nc); % remove padding bits' LLRs
    
    LLR_tikh   =   LLR_tikh';
    LLR_tikh   =   - LLR_tikh(:);%% negative sign to accomodate matlab
    LLR_tikh   =   LLR_tikh( 1: Nc); % remove padding bits' LLRs
    
    % Decoder
    dec_out_LLR   = step(hDec, LLR);
    dec_out_LLR   = dec_out_LLR(1:Nb); %% remove parity bits
    receivedBits   = ( dec_out_LLR<0 );
    
    dec_out_LLR_tikh   = step(hDec, LLR_tikh);
    dec_out_LLR_tikh   = dec_out_LLR_tikh(1:Nb); %% remove parity bits
    receivedBits_tikh   = ( dec_out_LLR_tikh<0 );
    % Error Count
    difference     = ( receivedBits ~= data);
    count        = count + sum( difference )/Nb;
    PER        = PER + any( difference );
    
    difference2     = ( receivedBits_tikh ~= data);
    count_tkh        = count_tkh + sum( difference2 )/Nb;
    PER        = PER + any( difference2 );
    %STAT(jj)   = sum( difference );  %%% number of bits in error
    %LLR_fliped = dec_out_LLR .* ( 2*data-1 ); %%
    %GMI(jj)    = 1- 1/Nb*sum( LLR2GMI( LLR_fliped  ) , 1);
end
BER(1)=count/BlocksN;
BER(2)=count_tkh/BlocksN;
PER=PER/BlocksN;
%1;