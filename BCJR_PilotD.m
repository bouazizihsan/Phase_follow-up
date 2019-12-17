function [BCJR]=BCJR_PilotD( Y, N0,sigw2, X,Lframe,payload,Pilot_Gamma_Apprx)
ss=size(Y);
Ns=ss(1);
D =ss(2);
Pilot_mean=unwrap( angle(Y./X) );
for n=1:Ns
    switch Pilot_Gamma_Apprx
        case 'SP'
            gamma{1,n} = Pilot_mean(n,:);
            gamma{2,n} = diag(N0./(2*abs(X(n,:)).*abs(Y(n,:))));
            
        case 'LTmax'
            gamma{1,n} = Pilot_mean(n,:);
            gamma{2,n} = diag(N0./(2*abs(X(n,:)).^2));
    end
end
%% A L P H A alpha(:,n) <-  alpha_{n-1}(p_n)
alpha{1,1}=[1 1] ;
alpha{2,1}=diag(ones(1,D).*Inf);
for n=2:Ns
    [Mnl,S]= mult_Gauss( alpha{1,n-1},alpha{2,n-1}, gamma{1,n-1}, gamma{2,n-1} );
    alpha{1,n}=Mnl;
    alpha{2,n}=S+sigw2;
end

%% B E T A beta(:,n) <- beta_{n+1}(p_n)
beta{1,Ns}=[1 1] ;
beta{2,Ns}=diag(ones(1,D).*Inf);
for n=Ns-1:-1:1
    [Mnl,S]= mult_Gauss( beta{1,n+1},beta{2,n+1}, gamma{1,n+1}, gamma{2,n+1} );
    beta{1,n}=Mnl;
    beta{2,n}=S+sigw2;
end

%% Combine
appos=[];
for n=1:Ns
[TT_ext{1,n},TT_ext{2,n}] = mult_Gauss( alpha{1,n},alpha{2,n}, beta{1,n},beta{2,n} );
[TT_apo{1,n},TT_apo{2,n}] = mult_Gauss( TT_ext{1,n},TT_ext{2,n}, gamma{1,n}, gamma{2,n});
appos=[appos;TT_apo{1,n}];
[alpha_gamma{1,n},alpha_gamma{2,n}] = mult_Gauss( alpha{1,n},alpha{2,n},gamma{1,n}, gamma{2,n} );
[beta_gamma{1,n},beta_gamma{2,n}] = mult_Gauss( beta{1,n},beta{2,n} ,gamma{1,n}, gamma{2,n} );

end
BCJR.alpha = alpha;
BCJR.beta  = beta;
BCJR.gamma = gamma; % direct reading of the phase
BCJR.ext = TT_ext;  %% extrinsic (to get the estimate excluding the pilot)
BCJR.apo = appos; %TT_apo;  %% aposteriori (to get the estimate for the pilots)
BCJR.parameters2=beta_gamma;
BCJR.parameters1=alpha_gamma;

end
function [m,v]=mult_Gauss(Mnl,Snl ,m_alpha,v_alpha)
        
        ss=diag(Snl);
        if(ss==Inf)
            m= m_alpha;
            v= v_alpha;
        elseif(diag(v_alpha)==Inf)
            m= Mnl;
            v= Snl;
        else
            S =v_alpha + Snl;
            m =(Snl*inv(S)*m_alpha.' + v_alpha*inv(S)*Mnl.').';
            v = inv(inv(v_alpha) + inv(Snl)) ;
        end
        
    end
