function [BER,PER,delta]= BErr_phase_noise_GT(hChan, hEnc, hDec,SNRdB)
BCJR1=[];
BCJR2=[];
Nblocks_stop=50;
Model=hChan.Model;
D=hChan.D;
M=hChan.M;
m= log2(M);
LLR_METHOD=hChan.LLR_METHOD;
SUM_METHOD= hChan.SUM_METHOD;
Lframe= hChan.Lframe;
Tech=hChan.Technique;
Phase_Estimation =hChan.PhaseEstimation;
sigw2= hChan.sigw2;
Npaquet = hChan.Npaquet;
% SNRdB   = hChan.SNRdB;    % SNR in dB
SNRlin  = 10.^(SNRdB/10); % Linear SNR
N0      = 1./SNRlin;
%%% coding info
[Nk,Nc]=size(hEnc.ParityCheckMatrix);
Nb=Nc-Nk; %% number of information bits
count=0;
PER=0;
delta=0;
for j=1:Npaquet
    %     fprintf('%d eme Paquet\n',j);
    Data=logical(randi([0 1],Nb,1));
    data_enc=step(hEnc,Data);
    switch Phase_Estimation
        case 'Yes'
            Signal=[];
            Xpilots=[];
            Xcons=[];
            Z=[];
            
            for i=1:D
                hMap=hChan.hMap(i);
                [K,mm]  = size(hMap.pntk0);
                [S,Zz,Xpilot,pilots,payload]=Generate_Data(data_enc,hMap,Nc,mm,m,Lframe);
                Signal=[Signal S];
                Xpilots=[Xpilots Xpilot];
                Xcons=[Xcons hChan.hMap(i).X];
                
                Z=[Z Zz*sqrt(N0/2)];
                
            end
            %             switch Model
            %                 case 'Tikhonov'
            %                     fMu=0;
            %                     fKappa=1/sigw2;
            %                     VonMises =vmrand(fMu, fKappa,size(Signal));
            %                     Phase_noise  =cumsum(VonMises);
            %                 case 'Gauss'
            w   = randn(length(Signal),D)+ 0;
            Phase_noise  = sqrt(sigw2).*cumsum(w);     % Wiener process
            %             end
            receivedSignal  = Signal.*exp(1i*Phase_noise)+Z;
            
            %             CCC={'r.','b.','g.'};
            switch Tech
                
                case 'HARQ_II'
                    switch Model
                        case 'Tikhonov'
                            [BCJR1]=Phase_track_Tikhonov_Hard1( receivedSignal(:,1), N0, sigw2, Xcons(:,1), Xpilots(:,1),pilots,Lframe,hChan.AlphaGammaProj);
                            LLR=LLR_Tikhonov1(M,receivedSignal(payload,1),hChan.hMap(1),N0,BCJR1,SUM_METHOD,payload);%                     plot((1:length(Phase_noise)),Phase_noise(:,1),'b');
                            %                     hold on;
                            %                      plot((1:length(Phase_noise)),unwrap(angle(BCJR1.Z_beta+BCJR1.Z_alpha)),'r');
                        case 'Gauss'
                            
                            %                     [BCJR2]    = Phase_track_BCJR_Gauss3( receivedSignal(:,1), N0, sigw2, hChan.hMap(1), pilots, Xpilots(:,1));
                            %                     yy=receivedSignal(payload,1).*exp(-1j*(BCJR2.Pmean(payload)+BCJR2.phase_shift(payload)  )); % + BCJR1.phase_shift(payload)
                            %                     LLR =LLR_PN(yy, 1, N0, BCJR2.var(payload), hChan.hMap(1), LLR_METHOD, SUM_METHOD); %, BCJR1.Taw(payload));
                            
                            [BCJR1]    = Phase_track_BCJR_Gauss(receivedSignal(:,1), N0, sigw2, Xcons(:,1), Xpilots(:,1),pilots,Lframe,hChan.Gamma_Apprx,hChan.Pilot_Gamma_Apprx,hChan.AlphaGammaProj,Phase_noise(:,1));
                            yy=receivedSignal(payload,1).*exp(-1j*(BCJR1.Pmean(payload)+ BCJR1.phase_shift(payload))); %
                            LLR =LLR_PN(yy, 1, N0, BCJR1.Pvar(payload), hChan.hMap(1), LLR_METHOD, SUM_METHOD);
                            
                            Diff=abs(Phase_noise(:,1)-(BCJR1.Pmean+ BCJR1.phase_shift)); %+ BCJR1.phase_shift
                            Diff(Diff>pi)=Diff(Diff>pi)-2*pi;
                            delta=delta+sum(Diff.^2 + BCJR1.Pvar);
                        case 'DP'
                            [BCJR1]=DP( receivedSignal(:,1), N0, sigw2, Xcons(:,1),  Xpilots(:,1),pilots);
                            LLR =LLR_Discretisation(BCJR1.p_apoX(payload,:),hChan.hMap(1),SUM_METHOD);
                            %                             yy=receivedSignal(payload,1).*exp(-1j*BCJR1.Pmean(payload)); %
                            %                             LLR =LLR_PN(yy, 1, N0, BCJR1.Pvar(payload), hChan.hMap(1), LLR_METHOD, SUM_METHOD);
                        case 'Interp'
                            [BCJR1]=Phase_track_interp(receivedSignal(:,1), N0, sigw2, Xcons(:,1), Xpilots(:,1),pilots,Lframe,hChan.Gamma_Apprx,hChan.Pilot_Gamma_Apprx);
                            yy=receivedSignal(payload,1).*exp(-1j*(BCJR1.Pmean(payload))); %
                            LLR =LLR_PN(yy, 1, N0, BCJR1.Pvar(payload), hChan.hMap(1), LLR_METHOD, SUM_METHOD);
                            
                    end
                    
                case 'HARQ_I'
                    [BCJR1]=Phase_track_Tikhonov_Hard1( receivedSignal(:,1), N0, sigw2, Xcons(:,1), Xpilots(:,1),pilots,Lframe );
                    LLR=LLR_Tikhonov1(M,receivedSignal(payload,1),hChan.hMap(1),N0,BCJR1,SUM_METHOD,payload);
                case 'HARQ_Somme_LLR'
                    LLR=zeros(length(payload),m);
                    for i=1:D
                        [BCJR1] = Phase_track_BCJR_Gauss_App1( receivedSignal(:,i), N0, sigw2, Xcons(:,i), Xpilots(:,i),pilots,Lframe,Phase_noise(:,i),hChan.hMap(i),Data,hDec,Nc,Nb,0);
                        yy=receivedSignal(payload,i).*exp(-1j*(BCJR1.Pmean(payload)+BCJR1.phase_shift(payload)));
                        LLR =LLR+LLR_PN(yy, 1, N0, BCJR1.Pvar(payload), hChan.hMap(i), LLR_METHOD, SUM_METHOD );
                        
                    end
                    
                case 'HARQ_Jointe'
                    %                                         for i=1:D
                    %                                             [BCJR1]    = Phase_track_BCJR_Gauss(receivedSignal(:,i), N0, sigw2, Xcons(:,i), Xpilots(:,i),pilots,Lframe,hChan.Gamma_Apprx,hChan.Pilot_Gamma_Apprx,hChan.AlphaGammaProj,Phase_noise(:,i));
                    %
                    %                                             figure(i)
                    %                                             plot((1:length(Phase_noise)),(BCJR1.Pmean+BCJR1.phase_shift),'m');
                    %                                         end
                    
                    [BCJR1]    =Phase_track_BCJR_Gauss_D( receivedSignal, N0, sigw2,Xcons, Xpilots,pilots,Lframe,payload,hChan.Gamma_Apprx,hChan.Pilot_Gamma_Apprx,hChan.AlphaGammaProj);
                    
                    %                                         figure(1)
                    %                                         hold on,
                    %                                         plot((1:length(Phase_noise)),(BCJR1.uapox(1:2:end)),'r');
                    %                                         plot((1:length(Phase_noise)),Phase_noise(:,1),'b');
                    %
                    %                                         figure(2)
                    %                                         hold on,
                    %                                         plot((1:length(Phase_noise)),(BCJR1.uapox(2:2:end)),'r');
                    %                                         plot((1:length(Phase_noise)),Phase_noise(:,2),'b');
                    %
                    %                     [BCJR2]    = Phase_track_BCJR_Gauss_D( receivedSignal, N0, sigw2,Xcons, Xpilots,pilots,Lframe,payload,hChan.Gamma_Apprx,hChan.Pilot_Gamma_Apprx,'OneSideCorr');% [hChan.hMap(1).Xsort hChan.hMap(2).Xsort]
                    %
                    %                                         figure(1)
                    %                                         hold on,
                    %                                         plot((1:length(Phase_noise)),(BCJR2.uapox(1:2:end)),'g');
                    %                                         hold off
                    %                                         figure(2)
                    %                                         hold on,
                    %                                         plot((1:length(Phase_noise)),(BCJR2.uapox(2:2:end)),'g');
                    %                                         hold off
                    
                    moy=[BCJR1.alpha_beta{:,1}].';
                    moy=[moy(1:2:end) moy(2:2:end)];
                    yy=receivedSignal(payload,:).*exp(-1j*moy);
                    LLR =LLR_PN_D(yy, diag([N0 N0]),BCJR1, hChan.hMap(1),hChan.hMap(2),payload,LLR_METHOD, SUM_METHOD );
                    %                     LLR=LLR_Discretisation(BCJR1.p_apo,hChan.hMap(1),hChan.hMap(2),SUM_METHOD );
                    
                    %                     LLR=zeros(length(payload),m);
                    %                     for i=1:D
                    %                         yy=receivedSignal(payload,i).*exp(-1j*(BCJR1.Pmean(:,i)));
                    %
                    %                         LLR =LLR+LLR_PN(yy, 1, N0, BCJR1.Pvar(:,i), hChan.hMap(i), LLR_METHOD, SUM_METHOD );
                    %                     end
                    
                    
                    %                     [BCJR1]=Phase_track_Tikhonov_D( receivedSignal, N0, sigw2, Xcons, Xpilots,pilots);
                    %                     LLR=LLR_Tikhonov1D(M,receivedSignal(payload,:),hChan.hMap(1),hChan.hMap(2),N0,BCJR1,SUM_METHOD);
                    
                otherwise
                    error('unknown Transmission Technique ');
                    
            end
            %             LLR2   =   LLR2';
            %             LLR2   =   LLR2(:);
            
            LLR   =   LLR';
            LLR   =  - LLR(:);
            LLR   = LLR(1:Nc);
            dec_llr_ph=step(hDec,LLR);
            
            dec_llr_ph=dec_llr_ph(1:Nb);
            Binest_ph =  (dec_llr_ph<0) ;
            difference     = ( Binest_ph ~= Data);
            count        = count + sum( difference );
            PER        = PER + any( difference );
            
        case 'None'
            Signal=[];
            Xpilots=[];
            Xcons=[];
            Z=[];
            for i=1:D
                hMap=hChan.hMap(i);
                [K,mm]  = size(hMap.pntk0);
                [S,Zz,Xpilot,pilots,payload]=Generate_Data(data_enc,hMap,Nc,mm,m,Lframe);
                
                Signal=[Signal S];
                Xpilots=[Xpilots Xpilot];
                Xcons=[Xcons hChan.hMap(i).Xsort];
                
                Z=[Z Zz.*sqrt(N0/2)];
                
            end
            
            Phase_noise=sqrt(sigw2)*randn(length(Signal),D);
            receivedSignal  =   Signal.*exp(1i*Phase_noise)+Z;
            
            switch Tech
                case 'HARQ_I'
                    LLR =LLR_MQAM_ph(M,receivedSignal(payload,1),Xcons(:,1).',N0,sigw2,LLR_METHOD,SUM_METHOD);
                    
                case 'HARQ_Somme_LLR'
                    LLR=zeros(length(payload),m);
                    for i=1:D
                        LLR =LLR+LLR_MQAM_ph(M,receivedSignal(payload,i),Xcons(:,i).',N0,sigw2,LLR_METHOD,SUM_METHOD);
                        
                    end
                otherwise
                    error('unknown Transmission Technique  ');
            end
            
            LLR = LLR';
            LLR = LLR(:);
            LLR = -LLR(1:Nc);
            
            dec_llr_ph=step(hDec,LLR);
            
            dec_llr_ph=dec_llr_ph(1:Nb);
            Binest_ph =  (dec_llr_ph<0) ;
            difference     = ( Binest_ph ~= Data);
            count = count+(sum(difference));
            PER        = PER + any( difference );
        otherwise
            error('unknown Phase Estimation ');
    end
    if PER>=Nblocks_stop, fprintf('early stop\n'), break, end
    clear  -regexp ^BCJR ^difference ^Binest_ph ^dec_llr_ph
end
switch Phase_Estimation
    case 'Yes'
        switch Tech
            case 'HARQ_II'
                delta=delta./(j.*length(Phase_noise));
                BER=(count./(j.*length(Data)));
                PER=(PER./j).^2;
            otherwise
                delta=delta./(Npaquet.*length(Phase_noise));
                BER=count./(Npaquet.*length(Data));
                PER=PER./Npaquet;
        end
        
    case 'None'
        switch Tech
            case 'HARQ_I'
                
                BER=(count./(Npaquet.*length(Data)));
                PER=(PER./Npaquet);
            case 'HARQ_II'
                
                BER=(count./(Npaquet.*length(Data)));
                PER=(PER./Npaquet);
            otherwise
                BER=count./(Npaquet.*length(Data));
                PER=PER./Npaquet;
        end
    otherwise
        error('unknown Phase Estimation ');
end
end