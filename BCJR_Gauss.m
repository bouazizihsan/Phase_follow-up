function [BCJR]=BCJR_Gauss( y, N0, sig2, X, pilot_symbols,pilot,payload)
ss           =size(y);
D            =ss(2);
Ns           =ss(1);
pilot_pos    = zeros(Ns,D);
pilot_pos(pilot)= pilot_symbols;
M            =length(X);
N0_2           =N0.*eye(D);
sigw2        =sig2.*eye(D);
Np           = 2^9;

Phase_limit  =pi;
dpp          = 2.*Phase_limit./Np;
pp           = linspace(-Phase_limit,Phase_limit,Np)';

pdf_Gauss=@(t,m,s2,d) real(exp(-0.5*sum(((t-m )'*inv(s2).*(t-m ).'),2) )./((2*pi).^(d/2)*det(s2).^(0.5)));

pdf_Gauss_log=@(t,m,s2,d) ((-0.5*(t-m )'*inv(s2)*(t-m )) - log( (2*pi).^(d/2)*det(s2).^(0.5) ));
g=pdf_Gauss( pp.' ,0 ,sigw2,D);
gamma=zeros(Np,Ns);
 gammaX=zeros(Np,Ns,M);

for n=1:Ns
    if( norm(pilot_pos(n))~=0)
        X_Hard=pilot_pos(n);
    else
        
        X_Hard=X;
        
    end
    MM=size(X_Hard,1);
    for i=1:MM  
        auxil= pdf_Gauss( (repmat(X_Hard(i),Np^D,1).*exp(1j.*pp)).' ,repmat(y(n,:),Np^D,1).', N0_2./2,D);
        gammaX(:,n,i)=auxil; 
        gamma(:,n)=gamma(:,n)+auxil;
    end
end


%ALPHA
alpha=zeros(Np,Ns);
alpha(:,:,1)=1;

for n=2:Ns
    
   
        alpha(:,n)= ifft( fft( alpha(:,n-1).* gamma(:,n-1), Np ) .* fft(fftshift(g),Np));
        alpha(:,n)=alpha(:,n)./(sum(sum(alpha(:,n))));
   
end

%Beta
beta=zeros(Np,Ns);
beta(:,Ns)=1;

for n=Ns-1:-1:1

        beta(:,n)= ifft( fft( beta(:,n+1).* gamma(:,n+1),Np ) .* fft(fftshift(g),Np)  );
        beta(:,n)=beta(:,n)./(sum(sum(beta(:,n))));

end


p_apo=repmat(beta.*alpha,1,1,M).*gammaX;
p_apo=p_apo./repmat(sum(sum(p_apo,1),3),Np,1,M);
p_apo=(sum(p_apo));

BCJR.p_apo=p_apo(:,payload,:);
end
