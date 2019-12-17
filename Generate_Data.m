
function [Signal,Z,Xpilots,pilots,payload]=Generate_Data(data_enc,hMap,Nc,mm,m,Lframe)
M=2.^m;
encodedDataPad  = [data_enc;randi([0 1],mm-rem(Nc,mm), 1)];% Padding bit
Ncodedsymbols=length(encodedDataPad)/mm;
WholeData= reshape(encodedDataPad,mm,Ncodedsymbols);
if mm~=m
    WholeData=[WholeData; randi([0 1], m-mm, Ncodedsymbols)]
end

encodeDataDec = bi2de( WholeData', 'left-msb' ); %% make it decimal
modSignal      = hMap.Xsort(encodeDataDec+1);
modSignal=modSignal.';
Ns=length(modSignal);
%% duration of the frame (payload+ one pilot)
Npilots=ceil(Ns/(Lframe-1))+1;
pilots  = (1:Lframe: Ns+Npilots )';
if(pilots(end)<Ns+Npilots)
    pilots=[pilots;Ns+Npilots];
end
payload=[1:Ns+Npilots];
payload(pilots)=[];
payload=payload';

%% duration of the frame (payload)

% Npilots=ceil((Ns)/(Lframe))+1;
% pilots  = (1:Lframe+1:(Npilots-1)*Lframe +Npilots)';
% if(pilots(end)~= Ns+Npilots)
%     pilots(end)=Ns+Npilots;
% end
% payload=[1:Ns+Npilots];
% payload(pilots)=[];
% payload=payload';

%% Generate the pilots
Signal=zeros(Npilots+Ns,1);%% Vect colonne
[hMapPilot] = create_constellation(2,'Gray','QAM'); 
var=randi([1 4], Npilots, 1);
Xpilots     = hMapPilot.Xsort(var); % random 4-QAM pilot symbols
Signal(payload)=modSignal;
Signal(pilots) = Xpilots; % symbols
Z=complex(randn(size(Signal)),randn(size(Signal)));

end