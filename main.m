%%main
% clear all;
% close all;
clc;
rng('shuffle');
SNRdB =[16:1:28]; % 64QAM
% SNRdB =[10:1:16 ];% 16QAM
%  SNRdB =3+[0:1:8];% 4QAM
% SNRdB =7+[0:1:8];% 4QAM sigw2 =0.05;
%  SNRdB =19+[0:1:12];% 64QAM sigw2 =0.05;
SAVE=1;
snrlin= 10.^(SNRdB./10); % 10 pour puissance et 20 log pour energie
hChan.Npaquet=1000;
% hChan.sigw2 =0.05;
hChan.sigw2 =0.0025;% 64QAM
% hChan.sigw2 =0.01;% 16QAM
%  hChan.sigw2 =0.01;% 4QAM

%%Model
hChan.Model = 'DP'; % Gauss, Tikhonov, DP, Interp
%% Independance:
%Envoie de meme symbole D fois et calculer les BERs
%indépendamment(Produit)
%% Somme_LLR :
%Envoie de meme symbole D fois et calculer les BERs en sommant
% Les LLRs(mémoriser les LLRs)

%% Jointe :
%Envoie de meme symbole D fois sur différants mapping et faire
%l'estimation de phase d'une façon jointe
Technique={'HARQ_II'}; %'HARQ_II' ,'HARQ_Somme_LLR', 'HARQ_Jointe'
%% modulation/demodulation
hChan.MOD='QAM'; %{ 'QAM' , 'PSK' };
Mapp= { {'Gray','Gray'}}; %,{'Gray','Ungerboeck'} {'Gray','Gray'},
hChan.M =64;  %%% constellation of size 2^m
% hChan.M =16;  %%% constellation of size 2^m
% hChan.M =4;  %%% constellation of size 2^m
hChan.D=2;
hChan.PhaseEstimation = 'Yes'; % None
AlphaGammaProj = {'Projection'}; % 'None' , 'Classtring' , 'Projection', 'OneSideCorr','BothSideCorr' ; 'None'= garder les différents mixtures


%% Encoding
% N=2016;
N=4032;
% N=16200;
%N=8064;
hChan.Lframe   = 17; % duration of the frame (payload+ one pilot)

switch N
    case 64800
        % encoder/decoder 64800
        Rate = 3/4; % 1/4, 1/3, 2/5, 1/2, 3/5, 2/3, 3/4, 4/5, 5/6, 8/9, and 9/10
        H  = dvbs2ldpc( Rate ); %parity matrix, code-length 64800
    case 16200
        % encoder/decoder 16200
        Rate = 8/9; % 1/5, 1/3, 2/5, 4/9, 3/5, 2/3, 11/15, 7/9, 37/45, 8/9
        H=t2_std_dvbt2_bl_ldpc( Rate, N, []);
    case {8064, 4032, 2016}
        Rate =3/4; % 7/8, 3/4 , 1/2
        setup.Rate = Rate;
        setup.Blength = N;
        setup.flag = 0; % 1
        
        Hobj = ldpc_hw_Hm(setup);
        H = Hobj.Hsprx;
end
hEnc    = comm.LDPCEncoder(H);
It      = 10; %[2 3 10]
%hDec    = comm.LDPCDecoder(H,'MaximumIterationCount', It,'DecisionMethod','Soft decision');
hDec    = comm.LDPCDecoder(H,'MaximumIterationCount', It,'DecisionMethod','Soft decision','OutputValue','Whole codeword'); %% get LLRs for all bits

CH_METHODS = {'SP1'} %CH_METHODS = {'LT','BLT','SP0','SP1','SP'};
Gamma_Apprx= {'LTmax'}; % ,'LT' 'SP',LTmax
hChan.Pilot_Gamma_Apprx='LTmax'; % 'SP'  'LTmax'
hChan.SUM_METHOD = 'sumexp';
for n=1:length(Technique)
    hChan.Technique=Technique{n};
    for nn=1:length(Gamma_Apprx)
        for kCH_Methods = 1:length(AlphaGammaProj)
            hChan.AlphaGammaProj=AlphaGammaProj{kCH_Methods}
            BER=zeros(1,length(SNRdB));
            PER=zeros(1,length(SNRdB));
            hChan.Gamma_Apprx=Gamma_Apprx{nn};
            hChan.hMap(1)= create_constellation(log2(hChan.M),Mapp{1}{1},hChan.MOD); % 'Gray', 'Ungerboeck',
            hChan.hMap(2) = create_constellation(log2(hChan.M),Mapp{1}{2},hChan.MOD);
            %     hChan.hMap(3) = create_constellation1(log2(hChan.M),Mapp{2}{2},hChan.MOD);
            hChan.LLR_METHOD = CH_METHODS{1};
            
            saving_file_name=['LDPC_',num2str(size(H,2))];
            saving_file_name=[saving_file_name,'_R=',num2str(Rate),'_It=',sprintf('%d_',It)];
            %%  M-QAM
            saving_file_name=[saving_file_name,'_',num2str(hChan.M),'_',hChan.MOD];
            %% Model
            saving_file_name=[saving_file_name,'_Model=', hChan.Model];
            switch hChan.Model
                case 'Gauss'
                    
                    %% LLR calculation method (for BICM)
                    saving_file_name=[saving_file_name,'_',hChan.hMap.MAPPING,'_',hChan.LLR_METHOD];
                    %% Gamma approximation method
                    saving_file_name=[saving_file_name,'_GammaApprx_',hChan.Gamma_Apprx];
                    %% Pilot Gamma approximation method
                    saving_file_name=[saving_file_name,'_PilotGammaApprx_',hChan.Pilot_Gamma_Apprx];
                    %% Pilot Gamma approximation method
                    saving_file_name=[saving_file_name,'_Correction_',hChan.AlphaGammaProj];
                    %% SNR range
                    saving_file_name=[saving_file_name,'_SNR=', num2str(SNRdB(1)),'dB_',num2str(SNRdB(end)),'dB'];
                    %% Npaquet
                    saving_file_name=[saving_file_name,'_Npaquet=',num2str(hChan.Npaquet)];
                    %% phase noise range
                    saving_file_name=[saving_file_name,'_sigw2=',num2str(hChan.sigw2)];
                    %% pilot¸
                    saving_file_name=[saving_file_name,'_L=',num2str(hChan.Lframe)];
                    %% phase estimation
                    saving_file_name=[saving_file_name,'_PhEst=', hChan.PhaseEstimation];
                    saving_file_name=[saving_file_name,'_XPilot_4QAM_'];
                case 'Tikhonov'
                    %% Pilot Gamma approximation method
                    saving_file_name=[saving_file_name,'_Correction=',hChan.AlphaGammaProj];
                    %% SNR range
                    saving_file_name=[saving_file_name,'_SNR=', num2str(SNRdB(1)),'dB_',num2str(SNRdB(end)),'dB'];
                    %% Npaquet
                    saving_file_name=[saving_file_name,'_Npaquet=',num2str(hChan.Npaquet)];
                    %% phase noise range
                    saving_file_name=[saving_file_name,'_sigw2=',num2str(hChan.sigw2)];
                    %% pilot¸
                    saving_file_name=[saving_file_name,'_L=',num2str(hChan.Lframe)];
                case 'DP'
                    %% LLR calculation method (for BICM)
                    saving_file_name=[saving_file_name,'_',hChan.hMap.MAPPING,'_',hChan.LLR_METHOD];
                    %% SNR range
                    saving_file_name=[saving_file_name,'_SNR=', num2str(SNRdB(1)),'dB_',num2str(SNRdB(end)),'dB'];
                    %% Npaquet
                    saving_file_name=[saving_file_name,'_Npaquet=',num2str(hChan.Npaquet)];
                    %% phase noise range
                    saving_file_name=[saving_file_name,'_sigw2=',num2str(hChan.sigw2)];
                    %% pilot¸
                    saving_file_name=[saving_file_name,'_L=',num2str(hChan.Lframe),'_2'];
                    
            end
            parfor i=1:length(SNRdB)
                fprintf('SNRdB = %d\n',SNRdB(i));
                [BER(i),PER(i),delta(i)]=BErr_phase_noise_GT(hChan, hEnc, hDec,SNRdB(i));
            end
            if SAVE
                parsave(['/Users/hsan/Documents/MATLAB/Res/64QAM/0.0025/'],saving_file_name,SNRdB,hChan,hEnc,hDec,BER,PER,delta)
            end
            clear -regexp ^BER ^PER ^delta ^saving_file_name3
        end
    end
    
end
disp('Gooooooooood. Thanks...')