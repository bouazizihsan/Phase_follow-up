%% example of the LDPC coded transmission

%------ D E F I N I T I O N S    of the general setup

%% Encoding
%N=64800;  % number of coded bits
%N=16200;
N=2016;
N=4032;
switch N
    case 64800
        % encoder/decoder 64800
        Rate = 8/9; % 1/4, 1/3, 2/5, 1/2, 3/5, 2/3, 3/4, 4/5, 5/6, 8/9, and 9/10
        H  = dvbs2ldpc( Rate ); %parity matrix, code-length 64800
    case 16200
        % encoder/decoder 16200
        Rate = 8/9; % 1/5, 1/3, 2/5, 4/9, 3/5, 2/3, 11/15, 7/9, 37/45, 8/9
        H=t2_std_dvbt2_bl_ldpc( Rate, N, []);
    case {8064, 4032, 2016}
        Rate = 1/2; % 7/8, 3/4 , 1/2
        setup.Rate = Rate;
        setup.Blength = N;
        setup.flag = 0; % 1
        
        Hobj = ldpc_hw_Hm(setup);
        H = Hobj.Hsprx;
end
hEnc    = comm.LDPCEncoder(H);
It      = 10; %[2 3 10]
%hDec    = comm.LDPCDecoder(H,'MaximumIterationCount', It,'DecisionMethod','Soft decision');
if length(It)>1
    for ii=1:length(It)
        hDec{ii}    = comm.LDPCDecoder(H,'MaximumIterationCount', It(ii),'DecisionMethod','Soft decision','OutputValue','Whole codeword'); %% get LLRs for all bits
    end
else
    hDec    = comm.LDPCDecoder(H,'MaximumIterationCount', It,'DecisionMethod','Soft decision','OutputValue','Whole codeword'); %% get LLRs for all bits
end

%%
hChan.Nblocks=1.e3 ; % maximum number of independent blocks (increase for low PER/BER)
hChan.Nblocks_stop=50;  % minimum number of blocks in error which trigger stop

% phase tracking
BCJR_PAR.GammaApproximation = 'SP'; %'SP', 'LTmax, % 'LTmax-F'%'LT0', %'BLT0',
BCJR_PAR.AlphaBetaPilotFilter = 'BothSidesCut'; % 'None', 'OneSideCut'  use only future/past pilots to remove the terms, 
                                              % 'BothSidesCut': use always future and past pilots to remove the terms
                                              % 'Multiply' : use future/past pilots to weight down the components

hChan.BCJR_PAR = BCJR_PAR;

% modulation/demodulation
hChan.M = 64;  %%% constellation of size M=2^m
hChan.hMap = create_constellation(log2(hChan.M),'Gray','QAM'); % 'Gray', 'Ungerboeck',
hChan.LLR_METHOD = 'SP'; %'LT'; %'BLT'  %'SP0'; %'SP1';

%% Channel and pilots
hChan.Lframe  = 17;  % 7, 9, 17, 49,%% payload + one pilot
hChan.Npilots = N/( (hChan.Lframe-1)*log2(hChan.M)) + 1;  %% number of pilots
%'P-Direct-Linear', 'P-BCJR-Linear', 'All-BCJR'; 'P-Perfect-Linear' ,
%'P-Direct-Linear-maxvar', 'Artificial-Parabola', 'All-Pilots', 'None'
hChan.PhaseEstimation = 'All-BCJR';

%------------------------------

% define the SNR range
%SNRdB=[5.6:0.2:6.2];  %% 4-QAM, rate 8/9, 64800
%SNRdB=15+[0:2:6];  %% 4-QAM, rate 8/9, 64800, sigw2=0.05
%SNRdB=7+[0:1:5];  %% 4-QAM, rate 7/8, 2016, sigw2=0.01, 'P-Direct-Linear'

%SNRdB=5+[0:1:5];  %% 4-QAM, rate 7/8, 2016, sigw2=0.01,
% SNRdB=5+[0:1:5];  %% 4-QAM, rate 7/8, 2016, sigw2=0.0025,

%SNRdB=10+[0:1:5];  %% 16-QAM, rate 3/4, 4032, sigw2=0.0025
%SNRdB=10+[0:1:5];  %% 16-QAM, rate 3/4, 4032, sigw2=0.01

%SNRdB=16+[0:1:6];  %% 64-QAM, rate 3/4, 4032, sigw2=0.0025
SNRdB=18+[0:1:8];  %% 64-QAM, rate 3/4, 4032, sigw2=0.01


%SNRdB = 32.4;% + [0:0.2:0.4];  %% 1024-QAM,
%SNRdB = 32 + [0:0.2:0.4];  %% 1024-QAM, rate 7/8, 8064, sigw2=4.e-5, P-Direct-Linear, pilot spacing L=17
%SNRdB = 35 + [0:0.2:0.6];  %% 1024-QAM, rate 7/8, 8064, sigw2=4.e-5, P-Direct-Linear, pilot spacing L=49
% define the range of the phase noise
%hChan.sigw2 = 3.78e-5; %[0.000];%,(0.001:0.002:0.03)];% used for 64QAM
%hChan.sigw2 = 0.01; %
%hChan.sigw2 = 0.0025; %

T0_sig2 = [0.05];
T1_GammaApproximation = {'SP'};%{'LT', 'SP'};
T2_LLR_METHOD = {'SP1'}; %{'LT', 'SP1'};
for t0=T0_sig2
for t1=1:length(T1_GammaApproximation)
for t2=1:length(T2_LLR_METHOD)
   hChan.sigw2 = t0;
   hChan.BCJR_PAR.GammaApproximation = T1_GammaApproximation{t1};
   hChan.LLR_METHOD = T2_LLR_METHOD{t2};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% saving file name for LDPC, rate,
saving_file_name=['./results/CODEC/','LDPC_',num2str(size(H,2))];
saving_file_name=[saving_file_name,'_R=',num2str(Rate),'_It=',sprintf('%d_',It)];
%%%  M-QAM
saving_file_name=[saving_file_name,'',num2str(hChan.M),'QAM'];
%%% LLR calculation method (for BICM)
saving_file_name=[saving_file_name,'_',hChan.hMap.MAPPING,'_',hChan.LLR_METHOD];
%%% SNR range
saving_file_name=[saving_file_name,'_SNR=', num2str(SNRdB(1)),'dB_',num2str(SNRdB(end)),'dB'];
%%% phase noise range
saving_file_name=[saving_file_name,'_sigmap2=',num2str(hChan.sigw2)];
%%% pilots
saving_file_name=[saving_file_name,'_L=',num2str(hChan.Lframe)];
%%% phase estimation
saving_file_name=[saving_file_name,'_PhEst=', hChan.PhaseEstimation];
%%% approximation of gamma
saving_file_name=[saving_file_name,'_GamApp=', hChan.BCJR_PAR.GammaApproximation];
%%% filter from pilots
saving_file_name=[saving_file_name,'_PltFilt=', hChan.BCJR_PAR.AlphaBetaPilotFilter];

BER  =  zeros( length(SNRdB), 1 );
PER  =  zeros( length(SNRdB), 1 );
GMI  = [];
PHI_stat = {};
%for gs=1:length(SIGMAp2)
gs=1;
fprintf(saving_file_name,'%s\n');
for gg=1:length(SNRdB)
    %hChan.SNRdB = SNRdB(gg);
    %%%% ---- everything happens here
    
    % performance under phase noise
    %[ BER(gg,gs), PER(gg,gs)] = codec_simulation_PN_run( hChan, hEnc, hDec );
    
    % tracking, phase=Wiener process
    [ BER(gg,gs), PER(gg,gs), GMI(:,gg)] = codec_simulation_phase_tracking_run( hChan, hEnc, hDec, SNRdB(gg) );
    
    % Turbo-tracking
    %[ BER(gg,gs), PER(gg,gs), GMI(:,gg)]=codec_simulation_turbo_phase_tracking_run( hChan, hEnc, hDec );
    
    % pure iid phase noise
    %     N0 = 10.^(-SNRdB(gg)/10);
    %     [ BER(gg,gs), PER(gg,gs)] = codec_simulation_run( N0, variance_PN, ...
    %                                      hEnc, hDec, hChan.hMap, hChan.LLR_METHOD );
    %%%%
    fprintf(';\n');
    
    %fprintf('%d/%d\n',gs,length(SIGMAp2));
    %%%  S A V E
end
save([saving_file_name,'.mat'],'SNRdB','BER','PER','GMI','hChan','hEnc','hDec')
end
end
end

