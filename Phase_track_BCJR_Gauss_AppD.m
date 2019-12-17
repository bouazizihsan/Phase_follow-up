function [BCJR]=Phase_track_BCJR_Gauss_AppD( y, N0, sig2, X, pilot_symbols,pilot,Lframe,payload,Phase_noise)
ss           =size(y);
D            =ss(2);
Ns           =ss(1);
pilot_pos    = zeros(Ns,2);
pilot_pos(pilot,:)= pilot_symbols;
M            =length(X);
N0_2           =N0.*eye(D);
sigw2        =sig2.*eye(D);
Np           = 2^8;
Npilots      = length(pilot_symbols);
Phase_limit  =3*sqrt(sig2*Lframe/4);  % we suppose the phase varies in the interval (-Phase_limit, Phase_limit)

dpp         = 2.*pi./Np;
p           = linspace(-pi,pi,Np).';

[M1,M2]=meshgrid(p);
% M1=M1.';
% M2=M2.';
M11=M1(:);
M22=M2(:);
pp=[M11 M22];

pdf_Gauss=@(t,m,s2,d) real(exp(-0.5*sum(((t-m )'*inv(s2).*(t-m ).'),2) )./((2*pi).^(d/2)*det(s2).^(0.5)));

pdf_Gauss_log=@(t,m,s2,d) ((-0.5*(t-m )'*inv(s2)*(t-m )) - log( (2*pi).^(d/2)*det(s2).^(0.5) ));

g= (reshape(pdf_Gauss( pp.' ,0 ,sigw2,D),Np,Np));

%processing
phase=angle(y(pilot,:)./pilot_symbols);
% for n=1:length(phase)-1
%    
%     deph= phase(n+1,:)-phase(n,:);
%     deph=deph(:);
%     fi=find(deph > pi); deph(fi)= deph(fi)-2.*pi;
%     fi=find(deph < -pi); deph(fi)= deph(fi)+2.*pi;
%     phase(n+1,:)=phase(n,:)+deph.';
% end
% new_phases=phase; %repmat(phase(1,:),Npilots,1) +[0 0;cumsum(dephasage)];

new_phases=unwrap(phase);
newphases=zeros(size(new_phases));
for n=1:length(new_phases)
    
    if(n==1)
            sigman_1=Inf;
            sigman1= (N0./(2*abs(pilot_symbols(n+1,:)).*abs(y(pilot(n)+Lframe,:))))+ Lframe.*sig2;
            sigman = (N0./(2*abs(pilot_symbols(n,:)).*abs(y(pilot(n),:)))) ;
            newphases(n,:) =(new_phases(n,:)./sigman + new_phases(n+1,:)./sigman1)./(1./sigman+1./sigman_1+1./sigman1);
        
        elseif(n==length(pilot))
            sigman1=Inf;
            sigman_1=(N0./(2*abs(pilot_symbols(n-1,:)).*abs(y(pilot(n)-Lframe,:))))+ Lframe.*sig2;
            sigman = (N0./(2*abs(pilot_symbols(n,:)).*abs(y(pilot(n),:)))) ;
            newphases(n,:) =(new_phases(n,:)./sigman + new_phases(n-1,:)./sigman_1)./(1./sigman+1./sigman_1+1./sigman1);
        
        else
            sigman_1=(N0./(2*abs(pilot_symbols(n-1,:)).*abs(y(pilot(n)-Lframe,:))))+ Lframe.*sig2;
            sigman1= (N0./(2*abs(pilot_symbols(n+1,:)).*abs(y(pilot(n)+Lframe,:))))+ Lframe.*sig2;
            sigman = (N0./(2*abs(pilot_symbols(n,:)).*abs(y(pilot(n),:)))) ;
        newphases(n,:)=(new_phases(n,:)./sigman + new_phases(n-1,:)./sigman_1 + new_phases(n+1,:)./sigman1)./(1./sigman+1./sigman_1+1./sigman1);
        
    end
        
    
end

phase_ref=interp1(pilot', newphases, (1:Ns)) ;

Theta_zero=zeros(M,D,Ns);
gamma=zeros(Np,Np,Ns);
for n=1:Ns
    auxi=[];
    aux1=[];
    aux_SNL=[];
    z=conj(y(n,:)).*exp(1j*Phase_noise(n,:));
    if( norm(pilot_pos(n,:))~=0)
        X_Hard=pilot_pos(n,:);
    else
        
        X_Hard=X;
        
    end  
    MM=size(X_Hard,1);
       Moy{1,n} =phase_ref(n,:)-imag(z.*(X_Hard))./abs(X_Hard).^2;
        for i=1:MM
        Snl_Xci{1,n} = N0/2.*abs(X_Hard).^2;
        Snl{1,n} = N0./(2*abs(X_Hard).^2);
        Xci_p0 {1,n}= pdf_Gauss( real( z*X_Hard ) ,abs(X_Hard).^2, Snl_Xci{1,n});
        % |X|*|Y|
        
        aux1=[aux1;diag(N0./(2*abs(X_Hard(i,:)).^2) )];
        aux_SNL=[aux_SNL;diag(N0./(2*abs(X_Hard(i,:)).*abs(y(n,:))))] ;
        auxi=[auxi;(pdf_Gauss_log( real(z.*(X_Hard(i,:))).',(abs(X_Hard).^2).',aux_SNL((i-1)*D+1:i*D,:),D) )];
        
        % |X|^2
        %         aux1=[aux1;diag(N0/2.*abs(X_Hard(i,:)).^2)];
        %         aux_SNL=[aux_SNL;diag(N0./(2*abs(X_Hard(i,:)).^2))] ;
        %         auxi=[auxi;(pdf_Gauss_log( real(y(n,:).*(X_Hard(i,:))).',(abs(X_Hard(i,:)).^2).',aux1((i-1)*D+1:i*D,:),D)')];

        
        end

    Theta_zero(:,:,n)=[Moy{1,n};zeros(M-size(Moy{1,n},1),ss(2))];
    Xci_p0{1,n}=auxi;
    Snl{1,n}=aux_SNL;

%     figure(12)
%     contour(M1,M2,gamma(:,:,n));
%     hold off
end

%Beta
beta=zeros(Np,Np,Ns);
beta(:,:,Ns)=1;
beta_p0{Ns,1}=[1 1] ;
beta_p0{Ns,2}=diag(ones(1,D).*Inf);

for n=Ns-1:-1:1
    if(n==241)
        stp=1;
    end
    [Mnl,S]=Projection_of_productD(Moy{1,n+1},Snl{1,n+1},Xci_p0{1,n+1},beta_p0{n+1,1},beta_p0{n+1,2},sigw2,Phase_limit,D,M1,M2,Phase_noise(n,:),dpp);
    
    beta_p0{n,1}= Mnl;
    beta_p0{n,2}=S;
end

%ALPHA
alpha=zeros(Np,Np,Ns);
alpha(:,:,1)=1;
alphafft=zeros(Np,Np,Ns);
alphafft(:,:,1)=1;
alpha_p0{1,1}=[1 1] ;
alpha_p0{1,2}=diag(ones(1,D).*Inf);
% alpha_p0{1,3}=[1 1];

for n=2:Ns
    if(n==96)
        stp=1;
    end
    [Mnl,S]=Projection_of_productD(Moy{1,n-1},Snl{1,n-1},Xci_p0{1,n-1},alpha_p0{n-1,1},alpha_p0{n-1,2},sigw2,Phase_limit,D,M1,M2,Phase_noise(n,:),dpp);
    
    alpha_p0{n,1}= Mnl;
    alpha_p0{n,2}=S;
    %     if(norm(pilot_pos(n,:))~=0)
    %         alpha_p0{n,2}=Snl{1,n};
    %     end
    
     alpha(:,:,n)= ifft2( fft2( alpha(:,:,n-1).* gamma(:,:,n-1), Np,Np ) .* fft2(fftshift(g), Np,Np)  );
%     alpha(:,:,n)= conv2( (alpha(:,:,n-1).* gamma(:,:,n-1)),g,'same');
    alpha(:,:,n)=alpha(:,:,n)./(sum(sum(alpha(:,:,n))));
%     figure(9)
%     contour(M1,M2,(alpha(:,:,n)));
%     hold on
%     plot(Phase_noise(n,1),Phase_noise(n,2),'r*');
%     hold off
%     
%     figure(10)
%     contour(M1,M2,gamma(:,:,n-1).*alpha(:,:,n-1));
%     hold on
%     plot(Phase_noise(n,1),Phase_noise(n,2),'r*');
%     hold off
end


p_apo=beta.*alpha;
p_apo=p_apo./repmat(sum(sum(p_apo,1)),Np,Np)./dpp;
for i=1:Ns
    Pmean(i,2)=angle(exp(1j.*p).'*sum(p_apo(:,:,i),2).*dpp);
    Pmean(i,1)=angle(exp(1j.*p).'*sum(p_apo(:,:,i),1)'.*dpp);
end
%Product of Alpha and Beta ==> one gaussian

M_alpha_beta=zeros(Ns,D);
V_alpha_beta=zeros(Ns*D,D);

for i=1:Ns
    
    if(i==1010)
        stp=1;
    end
    if(i==Ns)
        [Mnl,S]=Projection_of_productD(alpha_p0{i,1},alpha_p0{i,2},[0 0],beta_p0{i,1},beta_p0{i,2},0,Phase_limit,D,M1,M2,Phase_noise(i,:),dpp);
        
    else % cas ou il y a pas de varience inf
        [Mnl,S]=Projection_of_productD(beta_p0{i,1},beta_p0{i,2},[0 0],alpha_p0{i,1},alpha_p0{i,2},0,Phase_limit,D,M1,M2,Phase_noise(i,:),dpp);
        
        
    end
    
    V_alpha_beta((i-1)*D+1:i*D,:)= S;
    M_alpha_beta(i,:)= Mnl;
    
end



% Pmean=Pmean(payload,:);



p_apo=beta.*gamma.*alpha;
p_apo=p_apo./repmat(sum(sum(p_apo,1)),Np,Np)./dpp;
%  Pmean_apos1=M_alpha_beta(:,1); % alpha et beta c tout
%
%  Pmean_apos1=Pmean_apos(:,1);

q=[payload*D-1 payload*D].';
q=q(:);
V_alpha_beta=V_alpha_beta(q,:);

BCJR.Mnl=Theta_zero(:,:,payload);
BCJR.Pmean=M_alpha_beta(payload,:);
BCJR.Pvar=V_alpha_beta;

BCJR.Pmean_apos1=M_alpha_beta(:,1);

BCJR.Pmean_apos2=M_alpha_beta(:,2);
BCJR.discretisation=Pmean;
end
