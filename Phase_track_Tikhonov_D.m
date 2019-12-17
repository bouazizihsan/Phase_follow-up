function [BCJR]=Phase_track_Tikhonov_D( y, N0, sigw2, X, pilot_symbols,pilot,Phase_noise)

pdf_Gauss=@(t,m,s2,d) real(exp(-0.5*sum(((t-m )'*inv(s2).*(t-m ).'),2) )./((2*pi).^(d/2)*det(s2).^(0.5)));
D=2;
Ns           =length(y);
pilot_pos(pilot,:)= pilot_symbols;
Npilots      = length(pilot_symbols);
N0_2           =N0.*eye(D);
Np=2^8;
dpp         = 2.*pi./Np;
p           = linspace(-pi,pi,Np).';
[M1,M2]=meshgrid(p);
% M1=M1.';
% M2=M2.';
M11=M1(:);
M22=M2(:);
Npp=size(M11,1);
pp=[M11 M22];
g= (reshape(pdf_Gauss( pp.' ,0 ,sigw2,D),Np,Np));

gamma=zeros(Np,Np,Ns);
for n=1:Ns
   
    if( norm(pilot_pos(n,:))~=0)
        X_Hard=pilot_pos(n,:);
    else
        
        X_Hard=X;
        
    end
    MM=size(X_Hard,1);
    for i=1:MM
        
         gamma(:,:,n)=gamma(:,:,n) + reshape(pdf_Gauss( (repmat(X_Hard(i,:),Np^D,1).*exp(1j.*pp)).' ,repmat(y(n,:),Np^D,1).', N0_2./2,D),Np,Np);
        
    end
    gamma(:,:,n)=gamma(:,:,n)./(sum(sum(gamma(:,:,n))));
end

%ALPHA
alpha=zeros(Np,Np,Ns);
alpha(:,:,1)=1;
alpha_plot=zeros(Npp,1);
Z_alpha=zeros(Ns,2);
corr=1;
for n=2:Ns
    
    if(pilot_pos(n-1,:)~=0)
        
        X_Hard=pilot_pos(n-1,:);
    else
        X_Hard=X;
    end
    MM=size(X_Hard,1);
    if( (MM==1) && ((1-corr)>10^(-5)))
        z=[Z_alpha(n-1,:);[0 0]] + repmat((2/N0).*y(n-1,:).*conj(X_Hard),2,1);
         
        lambda= (abs(z(1,:))-abs(Z_alpha(n-1,:)))-0.5*log(abs(z)./repmat(abs(Z_alpha(n-1,:)),MM,1))-(abs(X_Hard).^2)./N0;
        lambda=[lambda; ( abs(z(2,:))-0.5*log(abs(z)./abs(Z_alpha(n-1,:)))-(abs(X_Hard).^2)./N0 )];
         lambda=[corr;1-corr].*exp(sum(lambda,2)-jac_logH((sum(lambda,2)).'));%% marginalisation
         corr=1;
    else
    z=repmat(Z_alpha(n-1,:),MM,1) + (2/N0).*repmat(y(n-1,:),MM,1).*conj(X_Hard);
     if(n==2)
        lambda=abs(z)-0.5*log(abs(z))-(abs(X_Hard).^2/N0);
    else
        lambda= (abs(z)-repmat(abs(Z_alpha(n-1,:)),MM,1))-0.5*log(abs(z)./repmat(abs(Z_alpha(n-1,:)),MM,1))-(abs(X_Hard).^2)./N0;
     end
     lambda=exp(sum(lambda,2)-jac_logH((sum(lambda,2)).'));%% marginalisation
  end
   
    
    lambda=exp(sum(lambda,2)-jac_logH((sum(lambda,2)).'));%% marginalisation
    % %     lambda= exp(lambda-max(lambda)-log( sum( exp(lambda-max(lambda))) )  );
    z_gamma =gamma_fct(sigw2,z);
    
    [zT,cor] = CMVM(lambda,z_gamma,D);
    corr=corr*cor;
    Z_alpha(n,:)=zT;
    
%     alpha_plot= exp(real(sum(repmat(Z_alpha(n,:),Npp,1).*exp(-1j*pp),2) ))./(2*pi*prod(besseli(0,abs(Z_alpha(n,:)) ),2) ) ;
%     alpha_plot=alpha_plot./sum(alpha_plot,1);
%     
%     alpha(:,:,n)= ifft2( fft2( alpha(:,:,n-1).* gamma(:,:,n-1), Np,Np ) .* fft2(fftshift(g), Np,Np)  );
%     alpha(:,:,n)=alpha(:,:,n)./(sum(sum(alpha(:,:,n))));
%     
%     figure(10)
%     contour(M1,M2,(alpha(:,:,n)));
%     hold on
%     plot(Phase_noise(n,1),Phase_noise(n,2),'r*');
%     hold off
%     alpha_plot=reshape(alpha_plot,[],length(M1));
%     figure(11)
%     contour(M1,M2,(alpha_plot));
%     hold on
%     plot(Phase_noise(n,1),Phase_noise(n,2),'r*');
%     hold off
end

%Beta
beta=zeros(Np,Np,Ns);
beta(:,:,Ns)=1;
beta_plot=zeros(Npp,1);
Z_beta=zeros(Ns,2);

for n=Ns-1:-1:1
    if(pilot_pos(n+1,:) ~= 0)
        X_Hard=pilot_pos(n+1,:);
    else
        X_Hard=X;
    end
    
    MM=size(X_Hard,1);
    z=repmat(Z_beta(n+1,:),MM,1) + (2/N0).*repmat(y(n+1,:),MM,1).*conj(X_Hard);
    
    if(n==(Ns-1))
        lambda=abs(z)-0.5*log(abs(z))-(abs(X_Hard).^2/N0);
    else
        lambda=(abs(z)-repmat(abs(Z_beta(n+1,:)),MM,1))-0.5*log(abs(z)./repmat(abs(Z_beta(n+1,:)),MM,1))-(abs(X_Hard).^2)/N0;
    end
    
    lambda=exp(sum(lambda,2)-jac_logH((sum(lambda,2)).'));
    %     lambda= exp(lambda-max(lambda)-log( sum( exp(lambda-max(lambda))) )  );
    z_gamma=gamma_fct(sigw2,z);
    
    [zT] = CMVM(lambda,z_gamma);
    Z_beta(n,:)=zT;
    
    beta_plot= exp(real(sum(repmat(Z_beta(n,:),Npp,1).*exp(-1j*pp),2) ))./(2*pi*prod(besseli(0,abs(Z_beta(n,:)) ),2) ) ;% + fact_inf./(2*pi); z=0 --> pdf=1./(2*pi)
    beta_plot=beta_plot./sum(beta_plot,1);
    
    beta(:,:,n)= ifft2( fft2( beta(:,:,n+1).* gamma(:,:,n+1), Np,Np ) .* fft2(fftshift(g), Np,Np)  );
    beta(:,:,n)=beta(:,:,n)./(sum(sum(beta(:,:,n))));
    if(n==145)
        stp=1;
    end
    figure(10)
    contour(M1,M2,(beta(:,:,n)));
    hold on
    plot(Phase_noise(n,1),Phase_noise(n,2),'r*');
    hold off
    beta_plot=reshape(beta_plot,[],length(M1));
    figure(11)
    contour(M1,M2,beta_plot);
    hold on
    plot(Phase_noise(n,1),Phase_noise(n,2),'r*');
    hold off
end

% BCJR.Z_betaP=Z_beta;
% BCJR.Z_alphaP=Z_alpha;
Z_alpha_beta=Z_beta + Z_alpha;
% Z_beta(pilot,:)=[];
% Z_alpha(pilot,:)=[];
Z_alpha_beta(pilot,:)=[];
% BCJR.Z_beta=Z_beta;
% BCJR.Z_alpha=Z_alpha;
BCJR.Z_alpha_beta=Z_alpha_beta;

end
function z=gamma_fct(sig,z)

z=z./(1+abs(z).*sig);

end
