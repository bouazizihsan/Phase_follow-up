
clc;
warning off;
snrdb = [0:2:20];
snrLin = 10.^(snrdb/10);
sigw2=0.01;
M=4;
mm=log2(M);
Npaquet=1;
BCJR=[];
BCJR2=[];

CH_METHODS = {'SP'};
%CH_METHODS = {'LT','BLT','SP0','SP1','SP','SPit','EX'};
SUM='sumexp';
Lframe   = 20; % duration of the frame (payload+ one pilot)

N=2016;
%N=16200;
%N=8064;

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

[Nk,Nc]=size(hEnc.ParityCheckMatrix);
Nb=Nc-Nk;
Ns=Nc/mm;
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
BER=zeros(length(CH_METHODS)*2,length(snrdb));
% BER=zeros(2,length(snrdb));
for kCH_Methods = 1:length(CH_METHODS)
LLR_METHOD = CH_METHODS{kCH_Methods}
for s= 1:length(snrdb)
    s
    % random data
    data  = logical(randi([0 1], Nb, 1));
    % Encode LDPC
    encodedData    = step(hEnc, data);
    WholeData = reshape( encodedData, mm, length(encodedData)/mm);
    sign= bi2de( WholeData.', 'left-msb' );
    
    [hMap]      = create_constellation(mm,'Gray'); % pilots constellation is 4-QAM
    XCons       = hMap.Xsort;
    XCons       = XCons.' ;
    x(payload)  = XCons(sign+1);
    
%     XCons=qammod(0:M-1,M,0,'gray'); %save constellation
%     power=sum(abs(XCons).^2)./length(XCons) ;
%     XCons= XCons./sqrt(power);
%     x(payload)= XCons(sign+1);
    
    % Generate the pilots
    Xpilots     = hMap.X(randi([1 M], Npilots, 1)); % random 4-QAM pilot symbols
    x(pilots) = Xpilots; % symbols

    w   =   sqrt(sigw2)*randn(size(x))+ 0; % sigma_p.^2: variance et 0= mean
    P  =   cumsum(w);     % Wiener process
    
    
    snr= snrLin(s);
    N0=1./(snr); % N0 total
    count=0;
    count_tikh=0;
    for h=1:Npaquet
        
        z=complex(randn(size(x)),randn(size(x)))./sqrt((2.*snr));% sigma
        y_ph=exp(1j*P).*x +z ;
 
        [BCJR]  = Phase_track_Tikhonov( y_ph, N0, sigw2, XCons, Xpilots,pilots,NClass);
        [BCJR2] = Phase_track_BCJR_Gauss_AppL( y_ph, N0, sigw2, XCons, Xpilots,pilots,NClass );
        
        LLR_tikh =LLR_Tikhonov(M,y_ph(payload),XCons,N0,BCJR.Z_alpha,BCJR.Z_beta,BCJR.Fact_alpha,BCJR.Fact_beta,SUM);
        LLR =LLR_cluster_Gauss(M,y_ph(payload),XCons,N0,BCJR2.L_alpha,BCJR2.L_beta,BCJR2.parameters1,BCJR2.parameters2,SUM,LLR_METHOD);
        %LLR =LLR_cluster_Gauss_Mat(M,y_ph(payload),XCons,N0,BCJR2.parameters1,BCJR2.parameters2,SUM,LLR_METHOD);
        ss=BCJR2.parameters1(:,2) + BCJR2.parameters2(:,2);
        var=BCJR2.parameters1(:,2).*BCJR2.parameters2(:,2)./(ss);
        u=(BCJR2.parameters1(:,2).*BCJR2.parameters2(:,1) + BCJR2.parameters1(:,1).*BCJR2.parameters2(:,2))./ss;
        yy=y_ph(payload).*exp(-1j*u);        
        if(s==8)
            abc=1;
        end
        LLR=LLR_PN(yy, 1, N0, var, hMap, LLR_METHOD, SUM );
                   
        LLR   =   LLR';
        LLR   =   - LLR(:);%% negative sign to accomodate matlab

        LLR_tikh   =   LLR_tikh';
        LLR_tikh   =   - LLR_tikh(:);%% negative sign to accomodate matlab
       
        % Decoder
        dec_out_LLR   = step(hDec, LLR);
        dec_out_LLR   = dec_out_LLR(1:Nb); %% remove parity bits
        receivedBits   = ( dec_out_LLR<0 );
        
        dec_out_LLR_tikh   = step(hDec, LLR_tikh);
        dec_out_LLR_tikh   = dec_out_LLR_tikh(1:Nb); %% remove parity bits
        receivedBits_tikh   = ( dec_out_LLR_tikh<0 );
        % Error Count
        difference     = ( receivedBits ~= data);
        count        = count + sum( difference );
        %PER        = PER + any( difference );
        
        difference2     = ( receivedBits_tikh ~= data);
        count_tikh        = count_tikh + sum( difference2 );
    end
    BER(1,s)=count./(Nb.*Npaquet); % BER of Gaussian approximation
    BER(2,s)=count_tikh./(Nb.*Npaquet);% BER of Tikhonov approximation
end
end
save BER_Tikhonov_Gauss_LDPC_Class.mat
