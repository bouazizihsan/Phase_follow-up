
clc;
warning off;
snrdb = [0:2:22];
snrLin = 10.^(snrdb/10);
sigw2=0.005;
M=16;
mm=log2(M);
Npaquet=1;
Nb=2016;
Ns=Nb/mm;
BCJR=[];
BCJR2=[];
NClass=3;
MOD={'QAM'};
CH_METHODS = {'SP'};
%CH_METHODS = {'LT','BLT','SP0','SP1','SP','SPit','EX'};
SUM='sumexp';
Lframe   = 20; % duration of the frame (payload+ one pilot)

Npilots=floor((Ns-1)/(Lframe-1))+1;
pilots  = (1:Lframe: Ns+Npilots-1 )';
if(pilots(end)~= Ns+Npilots)
    pilots(end)=Ns+Npilots;
end
% Npilots=length( pilots );
payload=[1:Ns+Npilots];
payload(pilots)=[];
payload=payload';

x=zeros(Ns+Npilots,1);
XCons=zeros(M,1);
BER=zeros(length(CH_METHODS)*3,length(snrdb));
% BER=zeros(2,length(snrdb));

for kCH_Methods = 1:length(CH_METHODS)
    LLR_METHOD = CH_METHODS{kCH_Methods}
    for s= length(snrdb):length(snrdb)
        s
        % random data
        data  = logical(randi([0 1], Nb, 1));
        WholeData = reshape( data, mm, Nb/mm);
        sign= bi2de( WholeData.', 'left-msb' );
        
        [hMap]      = create_constellation(mm,'Gray',MOD{1}); % pilots constellation is 4-QAM
        XCons       = hMap.Xsort;
        XCons       = XCons.' ;
        x(payload)  = XCons(sign+1);
        
        % Generate the pilots
        Xpilots     = hMap.X(randi([1 M], Npilots, 1)); % random 4-QAM pilot symbols
        x(pilots) = Xpilots; % symbols
        
        w   =  randn(size(x))+ 0;
        P  =sqrt(sigw2)*cumsum(w);     % Wiener process
        snr= snrLin(s);
        N0=1./(snr); 
        count=0;
        count_tikh=0;
        count_G=0;
        for h=1:Npaquet
            
            z=complex(randn(size(x)),randn(size(x)))./sqrt((2.*snr));% sigma
            y_ph=exp(1j*P).*x +z ;
            
            [BCJR]  = Phase_track_Tikhonov( y_ph, N0, sigw2, XCons, Xpilots,pilots,Lframe,NClass);
             [BCJR2] = Phase_track_BCJR_Gauss_AppL( y_ph, N0, sigw2, XCons, Xpilots,pilots,Lframe,NClass );
                                                     
            %[BCJR2]=Phase_track_BCJR_Gauss_App3( y_ph, N0, sigw2, hMap, pilots, Xpilots)
            
            yy=y_ph(payload).*exp(-1j* BCJR2.phase_shift(payload));
%             LLR_G =LLR_cluster_Gauss_Mat(M,yy,XCons,N0,BCJR.G1,BCJR.G2,SUM,LLR_METHOD);
            LLR =LLR_cluster_Gauss(M,yy,XCons,N0,BCJR2.L_alpha,BCJR2.L_beta,BCJR2.parameters1,BCJR2.parameters2,SUM,LLR_METHOD);
            LLR_tikh =LLR_Tikhonov(M,yy,XCons,N0,BCJR.Z_alpha,BCJR.Z_beta,BCJR.Fact_alpha,BCJR.Fact_beta,SUM);
            
            %LLR_G=LLR_Tikhonov(M,yy,XCons,N0,BCJR.Z_alphaG,BCJR.Z_betaG,BCJR.Fact_alphaG,BCJR.Fact_betaG,SUM);
           
            LLR   =   LLR';
            LLR   =    LLR(:);
               
% %             LLR_G   =   LLR_G';
% %             LLR_G   =    LLR_G(:);
            
            LLR_tikh   =   LLR_tikh';
            LLR_tikh   =    LLR_tikh(:);
            
            receivedBits   = ( LLR > 0 );
            
% %             receivedBits_G   = ( LLR_G > 0 );
            
            receivedBits_tikh   = ( LLR_tikh > 0 );
            % Error Count
            difference     = ( receivedBits ~= data);
            count        = count + sum( difference );
              
            
            difference2     = ( receivedBits_tikh ~= data);
            count_tikh        = count_tikh + sum( difference2 );
            
% %             difference3     = ( receivedBits_G ~= data);
% %             count_G        = count_G + sum( difference3 );    
             
        end
        BER(1,s)=count./(Nb.*Npaquet); % BER of Gaussian approximation
        BER(2,s)=count_tikh./(Nb.*Npaquet);% BER Tikh
% %         BER(3,s)=count_G./(Nb.*Npaquet);% BER Tikh_G
        
    end
end
figure (2)
semilogy(snrdb(1:length(BER(1,:))),BER(1,:), 'g')
hold on
semilogy(snrdb(1:length(BER(1,:))),BER(2,:), 'r')
% % semilogy(snrdb(1:length(BER(1,:))),BER(3,:), 'b')
title('BER')
xlabel('SNR-dB')
ylabel('BER')
legend('Gauss','Tikh-G')
%  save BER_Tikhonov_Gauss_Class001.mat