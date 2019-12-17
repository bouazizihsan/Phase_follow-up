%% Plot PER 
clf
clc
% %%64QAM
% %sigw2=0.012
% load('./results/G.vs.Tikh/64QAM/0.012/LDPC_4032_R=0.75_It=10__64_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_LT_PilotGammaApprx_SP_Correction=Projection_SNR=16dB_22dB_Npaquet=1000_sigw2=0.012_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat');
% semilogy(SNRdB(1:length(PER)),sqrt(PER), 'r-')
% load('./results/G.vs.Tikh/64QAM/0.012/LDPC_4032_R=0.75_It=10__64_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=Projection_SNR=16dB_22dB_Npaquet=1000_sigw2=0.012_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
% hold on,semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g-')
% load('./results/G.vs.Tikh/64QAM/0.012/LDPC_4032_R=0.75_It=10__64_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=OneSideCorr_SNR=16dB_22dB_Npaquet=1000_sigw2=0.012_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
% semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g--')
% load('./results/G.vs.Tikh/64QAM/0.012/LDPC_4032_R=0.75_It=10__64_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=BothSideCorr_SNR=16dB_22dB_Npaquet=1000_sigw2=0.012_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
% semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g-.')
% load('/Users/hsan/Documents/MATLAB/Master/results/G.vs.Tikh/64QAM/0.012/LDPC_4032_R=0.75_It=10__64_QAM_Model=Tikhonov_Correction=Projection_SNR=16dB_22dB_Hsan_Npaquet=1000_sigw2=0.012_L=17.mat')
% semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b-')
% load('./results/G.vs.Tikh/64QAM/0.012/LDPC_4032_R=0.75_It=10__64_QAM_Model=Tikhonov_Correction=OneSideCorr_SNR=16dB_22dB_Hsan_Npaquet=1000_sigw2=0.012_L=17.mat')
% semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b--')
% load('./results/G.vs.Tikh/64QAM/0.012/LDPC_4032_R=0.75_It=10__64_QAM_Model=Tikhonov_Correction=BothSideCorr_SNR=16dB_22dB_Hsan_Npaquet=1000_sigw2=0.012_L=17.mat')
% semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b-.')
% %sigw2=0.01
% load('./results/G.vs.Tikh/64QAM/0.01/LDPC_4032_R=0.75_It=10__64_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_LT_PilotGammaApprx_SP_Correction=Projection_SNR=16dB_22dB_Npaquet=1000_sigw2=0.01_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat');
% semilogy(SNRdB(1:length(PER)),sqrt(PER), 'r-')
% load('./results/G.vs.Tikh/64QAM/0.01/LDPC_4032_R=0.75_It=10__64_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=Projection_SNR=16dB_22dB_Npaquet=1000_sigw2=0.01_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
% hold on, semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g-')
% load('./results/G.vs.Tikh/64QAM/0.01/LDPC_4032_R=0.75_It=10__64_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=OneSideCorr_SNR=16dB_22dB_Npaquet=1000_sigw2=0.01_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
% semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g--')
% load('./results/G.vs.Tikh/64QAM/0.01/LDPC_4032_R=0.75_It=10__64_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=BothSideCorr_SNR=16dB_22dB_Npaquet=1000_sigw2=0.01_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
% semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g-.')
% load('./results/G.vs.Tikh/64QAM/0.01/LDPC_4032_R=0.75_It=10__64_QAM_Model=Tikhonov_Correction=Projection_SNR=16dB_22dB_Hsan_Npaquet=1000_sigw2=0.01_L=17.mat')
% semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b-')
% load('./results/G.vs.Tikh/64QAM/0.01/LDPC_4032_R=0.75_It=10__64_QAM_Model=Tikhonov_Correction=OneSideCorr_SNR=16dB_22dB_Hsan_Npaquet=1000_sigw2=0.01_L=17.mat')
% semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b--')
% load('./results/G.vs.Tikh/64QAM/0.01/LDPC_4032_R=0.75_It=10__64_QAM_Model=Tikhonov_Correction=BothSideCorr_SNR=16dB_22dB_Hsan_Npaquet=1000_sigw2=0.01_L=17.mat')
% semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b-.')
%sigw2=0.0025
load('./results/G.vs.Tikh/64QAM/0.0025/LDPC_4032_R=0.75_It=10__64_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_LT_PilotGammaApprx_SP_Correction=Projection_SNR=16dB_22dB_Npaquet=1000_sigw2=0.0025_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat');
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'r-')
load('./results/G.vs.Tikh/64QAM/0.0025/LDPC_4032_R=0.75_It=10__64_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=Projection_SNR=16dB_22dB_Npaquet=1000_sigw2=0.0025_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
hold on, semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g-')
load('./results/G.vs.Tikh/64QAM/0.0025/LDPC_4032_R=0.75_It=10__64_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=OneSideCorr_SNR=16dB_22dB_Npaquet=1000_sigw2=0.0025_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g--')
load('./results/G.vs.Tikh/64QAM/0.0025/LDPC_4032_R=0.75_It=10__64_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=BothSideCorr_SNR=16dB_22dB_Npaquet=1000_sigw2=0.0025_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g-.')
load('./results/G.vs.Tikh/64QAM/0.0025/LDPC_4032_R=0.75_It=10__64_QAM_Model=Tikhonov_Correction=Projection_SNR=16dB_22dB_Hsan_Npaquet=1000_sigw2=0.0025_L=17.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b-')
load('./results/G.vs.Tikh/64QAM/0.0025/LDPC_4032_R=0.75_It=10__64_QAM_Model=Tikhonov_Correction=OneSideCorr_SNR=16dB_22dB_Hsan_Npaquet=1000_sigw2=0.0025_L=17.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b--')
load('./results/G.vs.Tikh/64QAM/0.0025/LDPC_4032_R=0.75_It=10__64_QAM_Model=Tikhonov_Correction=BothSideCorr_SNR=16dB_22dB_Hsan_Npaquet=1000_sigw2=0.0025_L=17.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b-.')

legend('Gam=LT,GamP=SP,LLR=SP1,sigw2 =0.012','Gam=SP,GamP=SP,LLR=SP1,sigw2 =0.012','Gam=SP,GamP=SP,LLR=SP1,OneSideCorr,sigw2 =0.012','Gam=SP,GamP=SP,LLR=SP1,BothSideCorr,sigw2 =0.012','Tikhonov,sigw2 =0.012','Tikhonov,OneSideCorr,sigw2 =0.012','Tikhonov,BothSideCorr,sigw2 =0.012','Gam=LT,GamP=SP,LLR=SP1,sigw2 =0.01','Gam=SP,GamP=SP,LLR=SP1,sigw2 =0.01','Gam=SP,GamP=SP,LLR=SP1,OneSideCorr,sigw2 =0.01','Gam=SP,GamP=SP,LLR=SP1,BothSideCorr,sigw2 =0.01','Tikhonov,sigw2 =0.01','Tikhonov,OneSideCorr,sigw2 =0.01','Tikhonov,BothSideCorr,sigw2 =0.01','Gam=LT,GamP=SP,LLR=SP1,sigw2 =0.0025','Gam=SP,GamP=SP,LLR=SP1,sigw2 =0.0025','Gam=SP,GamP=SP,LLR=SP1,OneSideCorr,sigw2 =0.0025','Gam=SP,GamP=SP,LLR=SP1,BothSideCorr,sigw2 =0.0025','Tikhonov,sigw2 =0.0025','Tikhonov,OneSideCorr,sigw2 =0.0025','Tikhonov,BothSideCorr,sigw2 =0.0025')
title('64QAM. R=3/4. 1000Paquets');
ylabel('PER')
xlabel('SnrdB')

%4QAM
figure
%sigw2=0.012
load('./results/G.vs.Tikh/4QAM/0.012/LDPC_4032_R=0.875_It=10__4_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_LT_PilotGammaApprx_SP_Correction=Projection_SNR=5dB_10dB_Npaquet=1000_sigw2=0.012_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'r-')
load('./results/G.vs.Tikh/4QAM/0.012/LDPC_4032_R=0.875_It=10__4_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=Projection_SNR=5dB_10dB_Npaquet=1000_sigw2=0.012_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
hold on, semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g-')
load('./results/G.vs.Tikh/4QAM/0.012/LDPC_4032_R=0.875_It=10__4_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=OneSideCorr_SNR=5dB_10dB_Npaquet=1000_sigw2=0.012_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g--')
load('./results/G.vs.Tikh/4QAM/0.012/LDPC_4032_R=0.875_It=10__4_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=BothSideCorr_SNR=5dB_10dB_Npaquet=1000_sigw2=0.012_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g-.')
load('/Users/hsan/Documents/MATLAB/Master/results/G.vs.Tikh/4QAM/0.012/LDPC_4032_R=0.875_It=10__4_QAM_Model=Tikhonov_Correction=Projection_SNR=5dB_10dB_Hsan_Npaquet=1000_sigw2=0.012_L=17.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b-')
load('./results/G.vs.Tikh/4QAM/0.012/LDPC_4032_R=0.875_It=10__4_QAM_Model=Tikhonov_Correction=OneSideCorr_SNR=5dB_10dB_Hsan_Npaquet=1000_sigw2=0.012_L=17.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b--')
load('/Users/hsan/Documents/MATLAB/Master/results/G.vs.Tikh/4QAM/0.012/LDPC_4032_R=0.875_It=10__4_QAM_Model=Tikhonov_Correction=BothSideCorr_SNR=5dB_10dB_Hsan_Npaquet=1000_sigw2=0.012_L=17.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b-.')
%sigw2=0.01
load('./results/G.vs.Tikh/4QAM/0.01/LDPC_4032_R=0.875_It=10__4_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_LT_PilotGammaApprx_SP_Correction=Projection_SNR=5dB_10dB_Npaquet=1000_sigw2=0.01_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'r-')
load('./results/G.vs.Tikh/4QAM/0.01/LDPC_4032_R=0.875_It=10__4_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=Projection_SNR=5dB_10dB_Npaquet=1000_sigw2=0.01_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
hold on, semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g-')
load('./results/G.vs.Tikh/4QAM/0.01/LDPC_4032_R=0.875_It=10__4_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=OneSideCorr_SNR=5dB_10dB_Npaquet=1000_sigw2=0.01_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g--')
load('./results/G.vs.Tikh/4QAM/0.01/LDPC_4032_R=0.875_It=10__4_QAM_Model=Gauss_GrayGray_SP1_GammaApprx_SP_PilotGammaApprx_SP_Correction=BothSideCorr_SNR=5dB_10dB_Npaquet=1000_sigw2=0.01_L=17_PhEst=Yes_XPilot_4QAM_HARQ_II.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'g-.')
load('./results/G.vs.Tikh/4QAM/0.01/LDPC_4032_R=0.875_It=10__4_QAM_Model=Tikhonov_Correction=Projection_SNR=5dB_10dB_Hsan_Npaquet=1000_sigw2=0.01_L=17.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b-')
load('./results/G.vs.Tikh/4QAM/0.01/LDPC_4032_R=0.875_It=10__4_QAM_Model=Tikhonov_Correction=OneSideCorr_SNR=5dB_10dB_Hsan_Npaquet=1000_sigw2=0.01_L=17.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b--')
load('./results/G.vs.Tikh/4QAM/0.01/LDPC_4032_R=0.875_It=10__4_QAM_Model=Tikhonov_Correction=BothSideCorr_SNR=5dB_10dB_Hsan_Npaquet=1000_sigw2=0.01_L=17.mat')
semilogy(SNRdB(1:length(PER)),sqrt(PER), 'b-.')

legend('Gam=LT,GamP=SP,LLR=SP1,sigw2 =0.012','Gam=SP,GamP=SP,LLR=SP1,sigw2 =0.012','Gam=SP,GamP=SP,LLR=SP1,OneSideCorr,sigw2 =0.012','Gam=SP,GamP=SP,LLR=SP1,BothSideCorr,sigw2 =0.012','Tikhonov,sigw2 =0.012','Tikhonov,OneSideCorr,sigw2 =0.012','Tikhonov,BothSideCorr,sigw2 =0.012','Gam=LT,GamP=SP,LLR=SP1,sigw2 =0.01','Gam=SP,GamP=SP,LLR=SP1,sigw2 =0.01','Gam=SP,GamP=SP,LLR=SP1,OneSideCorr,sigw2 =0.01','Gam=SP,GamP=SP,LLR=SP1,BothSideCorr,sigw2 =0.01','Tikhonov,sigw2 =0.01','Tikhonov,OneSideCorr,sigw2 =0.01','Tikhonov,BothSideCorr,sigw2 =0.01','Gam=LT,GamP=SP,LLR=SP1,sigw2 =0.0025','Gam=SP,GamP=SP,LLR=SP1,sigw2 =0.0025','Gam=SP,GamP=SP,LLR=SP1,OneSideCorr,sigw2 =0.0025','Gam=SP,GamP=SP,LLR=SP1,BothSideCorr,sigw2 =0.0025','Tikhonov,sigw2 =0.0025','Tikhonov,OneSideCorr,sigw2 =0.0025','Tikhonov,BothSideCorr,sigw2 =0.0025')
title('4QAM. R=7/8. 1000Paquets');
ylabel('PER')
xlabel('SnrdB')
