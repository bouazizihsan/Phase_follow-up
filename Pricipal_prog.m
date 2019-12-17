% 4 QAM
%clear all;
clc;
SNRdB =[20];
SAVE=1;
snrlin= 10.^(SNRdB./10); % 10 pour puissance et 20 log pour energie
hChan.Npaquet=1000;
hChan.sigw2 =0.01;
Class=[1];

MOD= { 'QAM' }; %{ 'QAM' , 'PSK' };
% modulation/demodulation
hChan.M = 4;  %%% constellation of size 2^m

hChan.PhaseEstimation = 'Yes';


%% Encoding
N=2016;
%N=16200;
%N=8064;
hChan.Lframe   = 20; % duration of the frame (payload+ one pilot)

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
        Rate = 7/8; % 7/8, 3/4 , 1/2
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

BER=zeros(2*length(Class),length(SNRdB));
PER=zeros(2*length(Class),length(SNRdB));
%CH_METHODS = {'LT','BLT','SP0','SP1','SP'};
CH_METHODS = {'SP'};
hChan.SUM_METHOD = 'sumexp';
for n=1:length(Class)
for nn=1:length(MOD)
for kCH_Methods = 1:length(CH_METHODS)
    hChan.MOD=MOD{nn};
    hChan.hMap = create_constellation(log2(hChan.M),'Gray',hChan.MOD); % 'Gray', 'Ungerboeck',
    hChan.NClass=Class(n);
    hChan.LLR_METHOD = CH_METHODS{kCH_Methods};
    for i=1:length(SNRdB)
        fprintf('SNRdB = %d\n',SNRdB(i));
         hChan.SNRdB = SNRdB(i);
        [BER(2*n-1:2*n,i),PER(2*n-1:2*n,i)]=BErr_phase_noise(hChan, hEnc, hDec);
    end
    
    saving_file_name=['LDPC_',num2str(size(H,2))];
    saving_file_name=[saving_file_name,'_R=',num2str(Rate),'_It=',sprintf('%d_',It)];
    %%%  M-QAM
    saving_file_name=[saving_file_name,'_',num2str(hChan.M),'_',hChan.MOD];
    %%% LLR calculation method (for BICM)
    saving_file_name=[saving_file_name,'_',hChan.hMap.MAPPING,'_',hChan.LLR_METHOD];
    %%% SNR range
    saving_file_name=[saving_file_name,'_SNR=', num2str(SNRdB(1)),'dB_',num2str(SNRdB(end)),'dB'];
    %%% Npaquet
    saving_file_name=[saving_file_name,'_Npaquet=',num2str(hChan.Npaquet)];
    %%% phase noise range
    saving_file_name=[saving_file_name,'_sigmap2=',num2str(hChan.sigw2)];
    %%% pilots
    saving_file_name=[saving_file_name,'_L=',num2str(hChan.Lframe)];
    %%% phase estimation
    saving_file_name=[saving_file_name,'_PhEst=', hChan.PhaseEstimation];
    saving_file_name=[saving_file_name,'_Class=', num2str(hChan.NClass)];
    if SAVE
        save(['./results/CODEC/',saving_file_name,'.mat'],'SNRdB','BER','PER','hChan','hEnc','hDec')
    end
end
end
end
figure (1)
semilogy(SNRdB(1:length(BER(1,:))),BER(1,:), 'g')
hold on
semilogy(SNRdB(1:length(BER(1,:))),BER(2,:), 'r')
% semilogy(SNRdB(1:length(BER(1,:))),BER(3,:), 'b->')
% semilogy(SNRdB(1:length(BER(1,:))),BER(4,:), 'k->')
% semilogy(SNRdB(1:length(BER(1,:))),BER(5,:), 'c-*')
% semilogy(SNRdB(1:length(BER(1,:))),BER(6,:), 'y-*')
legend('Gauss1','Tikh1');%,'Tikh1','tikh2','Gauss3','Tikh3')
title('BER (64 QAM)')
xlabel('SNR-dB sigw2 =0.01')
ylabel('BER')
hold off