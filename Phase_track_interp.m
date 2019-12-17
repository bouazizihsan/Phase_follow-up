function [BCJR]=Phase_track_interp( y, N0, sigw2, X, pilot_symbols,pilot,Lframe,Gamma_Apprx,Pilot_Gamma_Apprx,AlphaGammaProj)
Ns           =length(y);
pilot_pos    = zeros(Ns,1);
pilot_pos(pilot)= pilot_symbols;
M            =length(X);

NP=length(pilot);
[BCJRP]=BCJR_Pilot( y(pilot), N0, Lframe*sigw2, pilot_symbols,Pilot_Gamma_Apprx);
nn = rem( (0:Ns-1)'  ,Lframe ) ;  %% position relative the the previous pilot
previous_pilot = ceil( (1:Ns)'/Lframe ); %% index of the previous pilot
previous_pilot =previous_pilot-ones(size(previous_pilot));
previous_pilot(end) = previous_pilot(end-1);
n=(0:Ns-1)';
delta=n-previous_pilot.*Lframe;
new_phases=BCJRP.apo(:,1); % unwrap(angle(y(pilot,:)./pilot_symbols)); %   Phase_noise(pilot); %
% phase_shift=(interp1(pilot', new_phases, (1:Ns))).' ; %=m
m=(1-delta./Lframe).*new_phases(previous_pilot+ones(size(previous_pilot)),:) + delta./Lframe.*new_phases(previous_pilot+2*ones(size(previous_pilot)),:);
sigma2=(1-delta./Lframe).*delta.*sigw2;
sigma2(pilot)=BCJRP.apo(:,2);
BCJR.Pmean=m;
BCJR.Pvar=sigma2;